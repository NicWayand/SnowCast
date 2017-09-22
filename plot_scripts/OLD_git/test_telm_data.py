import os
import time

csta = 'BNS' 
data_dir = '/home/nwayand/hydrology_staff_readonly/Centre/Marmot Creek/Telemetry/FortressMountain'
if os.name=='posix': # Linux
    # Hack to force datetimes to display in GMT/UTC (numpy 1.11.1 has fixed this but other dependent modules (pynio) can't handel numpy 1.11.1)
    os.environ['TZ'] = 'CST'
    time.tzset()

for c in [1]:
    import os
    import numpy as np
    import pandas as pd
    from astropy.io import ascii
    import datetime

    c_header = 4 # Header lines
    c_column_line = 1 # line where column names start
    c_delimiter = ','

    # Load in file
    cfile = csta + '_0015.dat' # Format for 15 min file
    dat = ascii.read(os.path.join(data_dir,cfile),header_start=c_column_line,data_start=c_header,delimiter=c_delimiter,exclude_names='N/A')
    datain = pd.DataFrame(dat.as_array())

    # Replace -9999 with nan (recomended by netcdf)
    datain.replace(-9999,np.NaN,inplace=True)

    # Make TIMESTAMP the index
    datain['TIMESTAMP'] = datain['TIMESTAMP'].astype('datetime64[ns]')
    datain = datain.set_index('TIMESTAMP')
    
    print datain.index

    
    # Telem files saved in CST, convert to MST
    datain.index = datain.index + datetime.timedelta(hours=CST_to_MST)
    # Adjust to local time zone (i.e. from UTC to MST, local_time_offset should = -7)
    #ds['datetime'] = pd.to_datetime(ds.datetime.values) + datetime.timedelta(hours=local_time_offset)
    
    print datain.index

    # Import header info 
    headerinfo = pd.read_csv(os.path.join(data_dir,cfile),nrows=2,skiprows=1)
    units = headerinfo.loc[0,:].tolist() # Grab first row of dataframe (units)
    units = units[1:] # Remove first value which is the units of the timestamp
