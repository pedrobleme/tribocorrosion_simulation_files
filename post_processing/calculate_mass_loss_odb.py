from odbAccess import *
import csv

# --- CONFIGURATION ---
odb_path = 'contatov1.odb'             # ODB file name
step_name = 'Step-2'                   # Step name
instance_name = 'PART-1-1'             # Instance name   
element_set_name = 'PART-1-1_SET-22'         # Optional: Specify a subset (e.g., 'ELSET_ALL')

# --- MANUAL VOLUME INPUT ---

single_element_volume = 1.0          # [mm^3] Average volume per element
# --------------------------------

def extract_continuous_mass_loss(odb_name, step_name, inst_name, elset=None):
    try:
        odb = openOdb(path=odb_name)
    except Exception as e:
        print("Error opening ODB: {}".format(e))
        return

    step = odb.steps[step_name]
    instance = odb.rootAssembly.instances[inst_name]
    
    report_file = 'continuous_mass_loss_report.csv'
    
    with open(report_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Frame', 'Time', 'Sum_of_SDV2', 'Total_Volume_Lost'])
        
        print("Processing frames for continuous damage extraction...")
        
        for i, frame in enumerate(step.frames):
            frame_time = frame.frameValue
            
            # 1. Safely extract SDV2
            try:
                sdv2_field = frame.fieldOutputs['SDV2']
            except KeyError:
                if i == 0:
                    # Frame 0 often lacks SDVs, write zeros and skip
                    writer.writerow([i, frame_time, 0.0, 0.0])
                    continue 
                else:
                    print("Error: 'SDV2' missing in frame {}.".format(i))
                    continue

            # 2. Subset to the specific contact/wear region
            if elset:
                try:
                    region = instance.elementSets[elset]
                    sdv2_field = sdv2_field.getSubset(region=region)
                except KeyError:
                    print("Error: Element set '{}' not found.".format(elset))
                    return
            
            # 3. Sum the damage fractions across all elements in the set
            sdv2_values = sdv2_field.values
            total_damage_fraction = 0.0
            
            for val in sdv2_values:
                total_damage_fraction += val.data
            
            # 4. Calculate total volume lost
            total_volume_lost = total_damage_fraction * single_element_volume
            
            writer.writerow([i, frame_time, total_damage_fraction, total_volume_lost])
            
            if i % 10 == 0:
                print("Processed Frame {} at Time {:.4f}".format(i, frame_time))

    odb.close()
    print("Report generated: {}".format(report_file))

if __name__ == '__main__':
    extract_continuous_mass_loss(odb_path, step_name, instance_name, element_set_name)