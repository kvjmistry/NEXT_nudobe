# Script to read in the nexus files and smear the hits along the track.
# Add some diffusion
# Re-bin the resultant track
import sys
import numpy  as np
import pandas as pd
from collections import Counter
import time

# USAGE:
# python3 SmearEventsGeneral.py <name of nexus input file name (remove .h5 extension)> <Scale Factor> <DL> <DT> <bin>
# e.g. python3 SmearEventsGeneral.py /Users/mistryk2/Packages/nexus/ATPC_0nuBB 1 0.408 1.440 1

# Record the start time
start_time = time.time()

# Load in the hits
print("Loading hits")
print("Filename: ", sys.argv[1]+".h5")
hits = pd.read_hdf(sys.argv[1]+".h5", 'MC/hits')
parts = pd.read_hdf(sys.argv[1]+".h5", 'MC/particles')
print("Finished loading hits")

# init the RNG
rng = np.random.default_rng()

# Diffusion values desired
DL = float(sys.argv[3]) # mm / sqrt(cm)
DT = float(sys.argv[4]) # mm / sqrt(cm)

# This is the scaling amount of diffusion
# scaling factor is in number of sigma of the diffusion values
# can set to zero for no diffusion, or 1 for default diffusion scaling of 1 sigma
diff_scaling = float(sys.argv[2])

binsize = int(sys.argv[5])

# Create the bins
# Configure here based on the detector geometry dimentions [mm]
# set detector active range with min, max and bin width in mm
xmin=-3000
xmax=3000
xbw=binsize

ymin=-3000
ymax=3000
ybw=binsize

zmin=0
zmax=6000
zbw=binsize

# This shifts the z pos of the events so 0 is at anode
# can set this to zero
z_shift = (zmax-zmin)/2.0

# Mean energy per e-. This splits up each G4 into E_hit/E_mean electrons
E_mean = 25e-6 # [eV]

print("Scaling Factor: ", diff_scaling)
print("DL: ", DL, "mm/sqrt(cm)")
print("DT: ", DT, "mm/sqrt(cm)")
print("Binsize is:",binsize," mm" )

# ----------- ----------- ----------- ----------- ----------- ----------- -----

# bins for x, y, z
xbins = np.arange(xmin, xmax+xbw, xbw)
ybins = np.arange(ymin, ymax+ybw, ybw)
zbins = np.arange(zmin, zmax+zbw, zbw)

# center bins for x, y, z
xbin_c = xbins[:-1] + xbw / 2
ybin_c = ybins[:-1] + ybw / 2
zbin_c = zbins[:-1] + zbw / 2


df_smear = []

# Define a function to smear the geant4 electrons uniformly between the steps
# Each electron is sampled uniformly half-way to the next hit and
# half-way towards the previous hit (based on step length)
# The ends of the track are sampled in the forward direction only 
def generate_random(row):
    r0 = np.array([row['x'], row['y'], row['z']])
    r1 = np.array([row['x'] - row['dx1']/2.0, row['y'] - row['dy1']/2.0, row['z'] - row['dz1']/2.0]) # backward delta
    r2 = np.array([row['x'] + row['dx2']/2.0, row['y'] + row['dy2']/2.0, row['z'] + row['dz2']/2.0]) # forward delta
    
    # Randomly pick to either move the electron forward or backward from the hit
    sampled_direction = np.random.choice([1, 2], p=[0.5, 0.5])

    if (sampled_direction == 1):
        random_number = rng.uniform(0, 1)
        new_r = r0+random_number*(r1 - r0)
    else:
        random_number = rng.uniform(0, 1)
        new_r = r0+random_number*(r2 - r0)

    x_smear, y_smear, z_smear = new_r[0], new_r[1], new_r[2]

    # Apply some diffusion to the electron too if the scaling is non-zero
    if (diff_scaling != 0):
        x = row['x'] # mm
        y = row['y'] # mm
        z = row['z'] # mm
        sigma_DL = diff_scaling*DL*np.sqrt(z/10.) # mm  
        sigma_DT = diff_scaling*DT*np.sqrt(z/10.) # mm 
    
        xy = np.array([x, y])
        cov_xy = np.array([[sigma_DT, 0], [0, sigma_DT]])
        
        xy_smear = rng.multivariate_normal(xy, cov_xy, 1)
        x_smear = xy_smear[0, 0]
        y_smear = xy_smear[0, 1]
        z_smear = rng.normal(z, sigma_DL)

    return pd.Series([x_smear, y_smear, z_smear], index=['x_smear', 'y_smear', 'z_smear'])

# Print the number of events:
print("Number of events to process: ", len(hits.event_id.unique()))
min_event_id = min( hits.event_id.unique())

# ---------------
# Main event Loop
# ---------------
for index, e in enumerate(hits.event_id.unique()):
    print("On Event:", e - min_event_id)

    # Record the end time
    end_time = time.time()

    # Calculate and print the runtime
    runtime = end_time - start_time
    print(f"Runtime: {runtime:.4f} seconds")

    # Select the event
    event = hits[hits.event_id == e]
    event_part = parts[parts.event_id == e]
    
    # Shift z-values so 0 is at the anode
    event.z = event.z+z_shift

    # Calc number of electrons in a hit
    event["n"] = round(event["energy"]/E_mean)

    # Loop over the particles and get the differences between steps ------
    particles = event.particle_id.unique()

    smear_df = []

    for idx, p in enumerate(particles):

        # Get hits for particle i in the event
        temp_part = event[event.particle_id == p]
        particle_name = event_part[event_part.particle_id == p].particle_name.iloc[0]

        nrows = len(temp_part)

        # This dataframe contains the difference in distance between hits
        diff_df = temp_part[['x', 'y', 'z']].diff()
        diff_df.iloc[0] = 0
        extra_row = pd.DataFrame({'x': [0], 'y': [0], 'z': [0]})
        diff_df = pd.concat([diff_df, extra_row])
        diff_df = diff_df.rename(columns={'x': 'dx', 'y': 'dy', 'z': 'dz'})

        # We dont want to smear over the gamma steps
        # Only their daughter electrons
        if (particle_name == "gamma"):
            diff_df["dx"] = 0*diff_df["dx"]
            diff_df["dy"] = 0*diff_df["dy"]
            diff_df["dz"] = 0*diff_df["dz"]
        
        dx1 = [] # backward delta from hit
        dy1 = []
        dz1 = []
        dx2 = [] # forward delta from hit
        dy2 = []
        dz2 = []
        index_arr = []

        # Get arrays for the distance dataframes
        for index in range(len(diff_df)-1):
            dx1.append(diff_df.iloc[index].dx)
            dy1.append(diff_df.iloc[index].dy)
            dz1.append(diff_df.iloc[index].dz)
            dx2.append(diff_df.iloc[index+1].dx)
            dy2.append(diff_df.iloc[index+1].dy)
            dz2.append(diff_df.iloc[index+1].dz)

        index_arr = diff_df.index.to_numpy()
        index_arr = index_arr[:-1]

        # Make a dataframe of the deltas
        data = {
            'dx1': dx1,
            'dx2': dx2,
            'dy1': dy1,
            'dy2': dy2,
            'dz1': dz1,
            'dz2': dz2,
        }

        new_df = pd.DataFrame(data, index=index_arr)
        smear_df.append(new_df)

    # Concatenate DataFrames along rows (axis=0)
    smear_df = pd.concat(smear_df)

    # Now merge to the main df
    event = pd.merge(event, smear_df, left_index=True, right_index=True, how='inner')

    # Create a new DataFrame with duplicated rows, so we can smear each electron by diffusion
    electrons = pd.DataFrame(np.repeat(event[["event_id",'x', 'y', 'z', 'dx1', 'dx2', 'dy1', 'dy2', 'dz1','dz2']].values, event['n'], axis=0), columns=["event_id",'x', 'y', 'z', 'dx1', 'dx2', 'dy1', 'dy2', 'dz1','dz2'])

    # Reset the index of the new DataFrame if needed
    electrons = electrons.reset_index(drop=True)

    # Now apply some smearing to each of the electrons
    # Apply the function to create new columns
    new_columns     = electrons.apply(generate_random, axis=1)
    electrons_smear = pd.concat([electrons, new_columns], axis=1)
    electrons_smear["energy"] = E_mean # MeV

    # We need to set this to make sure we keep the information about the unbinned positions in the weighting
    electrons_smear['x'] = electrons_smear['x_smear']
    electrons_smear['y'] = electrons_smear['y_smear']
    electrons_smear['z'] = electrons_smear['z_smear']

    # Now lets bin the data
    electrons_smear['x_smear'] = pd.cut(x=electrons_smear['x_smear'], bins=xbins,labels=xbin_c, include_lowest=True)
    electrons_smear['y_smear'] = pd.cut(x=electrons_smear['y_smear'], bins=ybins,labels=ybin_c, include_lowest=True)
    electrons_smear['z_smear'] = pd.cut(x=electrons_smear['z_smear'], bins=zbins,labels=zbin_c, include_lowest=True)

    # Loop over the rows in the dataframe and sum the energies of all electrons in a bin. 
    # Also change the bin center to use the mean x,y,z position
    x_mean_arr = []
    y_mean_arr = []
    z_mean_arr = []
    energy_mean_arr = []
    x_mean_arr_temp = np.array([])
    y_mean_arr_temp = np.array([])
    z_mean_arr_temp = np.array([])
    summed_energy = 0
    event_id = 0

    counter = 0

    # Sort so all the bin labels are next to one another
    electrons_smear = electrons_smear.sort_values(by=['x_smear', 'y_smear', 'z_smear'])

    # Loop over all bins and aggregate to get total energy in each bin and their
    # mean x,y,z position
    for index, row in electrons_smear.iterrows():

        # First row 
        if (counter == 0):
            temp_x = row["x_smear"]
            temp_y = row["y_smear"]
            temp_z = row["z_smear"]
            summed_energy +=row["energy"]
            event_id = row["event_id"]
            x_mean_arr_temp = np.append(x_mean_arr_temp, row["x"])
            y_mean_arr_temp = np.append(y_mean_arr_temp, row["y"])
            z_mean_arr_temp = np.append(z_mean_arr_temp, row["z"])
            counter+=1
            continue

        # Final row
        if index == electrons_smear.index[-1]:
            x_mean_arr_temp = np.append(x_mean_arr_temp, row["x"])
            y_mean_arr_temp = np.append(y_mean_arr_temp, row["y"])
            z_mean_arr_temp = np.append(z_mean_arr_temp, row["z"])
            summed_energy +=row["energy"]

            if (summed_energy != 0): 
                x_mean_arr = np.append(x_mean_arr,np.mean(x_mean_arr_temp))
                y_mean_arr = np.append(y_mean_arr,np.mean(y_mean_arr_temp))
                z_mean_arr = np.append(z_mean_arr,np.mean(z_mean_arr_temp))
                energy_mean_arr.append(summed_energy)


        # Same bin
        if (row["x_smear"] == temp_x and row["y_smear"] == temp_y and row["z_smear"] == temp_z):
            x_mean_arr_temp = np.append(x_mean_arr_temp, row["x"])
            y_mean_arr_temp = np.append(y_mean_arr_temp, row["y"])
            z_mean_arr_temp = np.append(z_mean_arr_temp, row["z"])
            summed_energy +=row["energy"]

        # Aggregate and store for next 
        else:
            if (summed_energy != 0): 
                x_mean_arr = np.append(x_mean_arr,np.mean(x_mean_arr_temp))
                y_mean_arr = np.append(y_mean_arr,np.mean(y_mean_arr_temp))
                z_mean_arr = np.append(z_mean_arr,np.mean(z_mean_arr_temp))
                energy_mean_arr.append(summed_energy)
            
            temp_x = row["x_smear"]
            temp_y = row["y_smear"]
            temp_z = row["z_smear"]
            summed_energy = 0
            x_mean_arr_temp = np.array([])
            y_mean_arr_temp = np.array([])
            z_mean_arr_temp = np.array([])
            
            
            x_mean_arr_temp = np.append(x_mean_arr_temp, row["x"])
            y_mean_arr_temp = np.append(y_mean_arr_temp, row["y"])
            z_mean_arr_temp = np.append(z_mean_arr_temp, row["z"])
            summed_energy +=row["energy"]
            

        counter+=1

    events = np.ones_like(energy_mean_arr)*event_id

    # Make the dataframe again
    electrons_smear = pd.DataFrame({  "event_id" : events, "x" : x_mean_arr,  "y" : y_mean_arr,  "z" : z_mean_arr,  "energy" : energy_mean_arr  }) 


    # File writing
    electrons_smear = electrons_smear.sort_values(by=['event_id', 'z', 'x', 'y'])

    electrons_smear['event_id'] = electrons_smear['event_id'].astype(int)
    electrons_smear['z'] = electrons_smear['z'].astype('float32')
    electrons_smear['x'] = electrons_smear['x'].astype('float32')
    electrons_smear['y'] = electrons_smear['y'].astype('float32')
    # df['energy'] = df['energy']*1e6
    electrons_smear['energy'] = electrons_smear['energy'].astype('float32')

    df_smear.append(electrons_smear)


df_smear_merge = pd.concat(df_smear, ignore_index=True)

print("Saving events to file: ", sys.argv[1]+"_"+str(binsize)+"mm_smear.h5")
with pd.HDFStore(sys.argv[1]+"_"+str(binsize)+"mm_smear.h5", mode='w', complevel=5, complib='zlib') as store:
    # Write each DataFrame to the file with a unique key
    store.put('MC/particles', parts, format='table')
    store.put('MC/hits', df_smear_merge, format='table')

# Record the end time
end_time = time.time()

# Calculate and print the runtime
runtime = end_time - start_time
print(f"Runtime: {runtime:.4f} seconds")