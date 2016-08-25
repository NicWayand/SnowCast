def get_recent_CRHO_data_from_telem(csta,data_dir,output_dir):
    import os
    import numpy as np
    import pandas as pd
    from astropy.io import ascii
    
    # Move to data dir
    os.chdir(data_dir) 

    c_header = 4 # Header lines
    c_column_line = 1 # line where column names start
    c_delimiter = ','

    # Load in file
    cfile = csta + '_0015.dat' # Format for 15 min file
    dat = ascii.read(cfile,header_start=c_column_line,data_start=c_header,delimiter=c_delimiter,exclude_names='N/A')
    datain = pd.DataFrame(dat.as_array())

    # Replace -9999 with nan (recomended by netcdf)
    datain.replace(-9999,np.NaN,inplace=True)

    # Make TIMESTAMP the index
    datain['TIMESTAMP'] = datain['TIMESTAMP'].astype('datetime64[ns]')
    datain = datain.set_index('TIMESTAMP')

    # Import header info 
    headerinfo = pd.read_csv(cfile,nrows=2,skiprows=1)
    units = headerinfo.loc[0,:].tolist() # Grab first row of dataframe (units)
    units = units[1:] # Remove first value which is the units of the timestamp

    return {'datain':datain, 'units':units}
