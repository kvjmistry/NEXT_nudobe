import re
import glob as glob
import subprocess
import pandas as pd

# Example grep output

# Run the grep command
command = 'grep -ri "Killed" ./log/'
result = subprocess.run(command, shell=True, capture_output=True, text=True)

# Get the output from the command
grep_output = result.stdout

# Regex to extract the required numbers
pattern = r'GN_(\d+)_(\d+)\.err'
matches = re.findall(pattern, grep_output)

models = []
nmes = []
binnings = []
jobids = []
pressures = []
submissionids = []

for line in matches:
    # Read the first 4 lines of the file
    file_path = f"log/GN_{line[0]}_{line[1]}.out"
    with open(file_path, 'r') as file:
        lines = [file.readline().strip() for _ in range(5)]

    # Combine the lines into a single string
    file_content = "\n".join(lines)

    # Extract the model, NME, and binning using regex
    model_match = re.search(r'Model is:\s*(\w+)', file_content).group(1) if re.search(r'Model is:\s*(\w+)', file_content) else None
    nme_match = re.search(r'NME is:\s*(\w+)', file_content).group(1) if re.search(r'NME is:\s*(\w+)', file_content) else None
    binning_match = re.search(r'BINNING is:\s*([\w_]+)', file_content).group(1) if re.search(r'BINNING is:\s*([\w_]+)', file_content) else None
    
    pressure_string = re.search(r'H5File is:\s*(.+)', file_content).group(1) if re.search(r'H5File is:\s*(.+)', file_content) else None
    pressure_match = re.search(r'/(\d+bar)/', pressure_string)
    pressure_match = pressure_match.group(1) if pressure_match else None
    models.append(model_match)
    nmes.append(nme_match)
    binnings.append(binning_match)
    jobids.append(line[1])
    pressures.append(pressure_match)
    submissionids.append(line[0])


df = pd.DataFrame( {"Model" : models, "NME" : nmes, "Bin": binnings, "pressure": pressures, "JOBID": jobids, "subid": submissionids }  )
df["JOBID"] = df["JOBID"].astype(int)
df["subid"] = df["subid"].astype(int)


def PrintFailed(df, model, nme, pressure, Bin):
    temp = df[ (df["Model"] == model) & (df["NME"] == nme) & (df["Bin"] == Bin) & (df["pressure"] == pressure) ]
    if (len(temp) == 0):
        return
    print(f"{model}, {nme}, {pressure}, {Bin}")
    for sub_id in temp.subid.unique():
        temp2 = temp[temp.subid == sub_id]
        print(sub_id)
        temp2 = temp2.sort_values(by=["JOBID"])
        comma_separated = ','.join(temp2['JOBID'].astype(str))
        print(comma_separated)
    print("\n")

PrintFailed(df, "mbb", "SM", "1bar","1mm_smear")
PrintFailed(df, "mbb", "SM", "1bar","2mm_smear")
PrintFailed(df, "mbb", "SM", "1bar","4mm_smear")
PrintFailed(df, "mbb", "SM", "1bar","10mm_smear")

PrintFailed(df, "mbb", "SM", "5bar","1mm_smear")
PrintFailed(df, "mbb", "SM", "5bar","2mm_smear")
PrintFailed(df, "mbb", "SM", "5bar","4mm_smear")
PrintFailed(df, "mbb", "SM", "5bar","10mm_smear")

PrintFailed(df, "mbb", "SM", "10bar","1mm_smear")
PrintFailed(df, "mbb", "SM", "10bar","2mm_smear")
PrintFailed(df, "mbb", "SM", "10bar","4mm_smear")
PrintFailed(df, "mbb", "SM", "10bar","10mm_smear")

PrintFailed(df, "mbb", "SM", "15bar","1mm_smear")
PrintFailed(df, "mbb", "SM", "15bar","2mm_smear")
PrintFailed(df, "mbb", "SM", "15bar","4mm_smear")
PrintFailed(df, "mbb", "SM", "15bar","10mm_smear")