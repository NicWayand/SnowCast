import numpy as np

# Save Figure
def save_figure(f,file_out,fig_res):


    f.savefig(file_out,bbox_inches='tight',dpi=fig_res)

# Function that makes two data sets common:
def make_common(ds_obs, ds_mod, dt_eval, dt_eval_hr, remove_missing=True, percent_nan_allowed=20):


    # dt_eval = 'H' # MS (month start), H (hour), W (week)
    # How do we handel missing values within an aggregation period?
    # We use a threshold (percent_nan_allowed) for % missing within a period,
    skipna = True # Skips nans when doing aggregation, then we remove those periods with less than percent_nan_allowed

    # Common variables
    obs_vars = ds_obs.data_vars
    mod_vars = ds_mod.data_vars
    # Find common variables
    com_vars = np.sort(list(set(obs_vars).intersection(mod_vars)))
    # Extract only common vars
    ds_obs = ds_obs[com_vars]
    ds_mod = ds_mod[com_vars]
    print("Common variables are:", com_vars)
    print("")

    # Common stations
    obs_sta = ds_obs['station'].values
    mod_sta = ds_mod['station'].values
    # Find common variables
    com_sta = np.sort(list(set(obs_sta).intersection(mod_sta)))
    print("Common stations are:", com_sta)
    print("")
    if not len(com_sta):
        print(obs_sta)
        print(mod_sta)
        raise ValueError('No common stations found')

    # Extract only common stations
    ds_obs = ds_obs.sel(station=com_sta)
    ds_mod = ds_mod.sel(station=com_sta)

    # Pre-trim time to common time start and end (speeds up time step aggregation)
    T_start = np.max([ds_mod.time.values[0], ds_obs.time.values[0]])
    T_end = np.min([ds_mod.time.values[-1], ds_obs.time.values[-1]])
    print(T_start, T_end)
    ds_obs = ds_obs.where((ds_obs.time >= T_start) & (ds_obs.time <= T_end), drop=True)
    ds_mod = ds_mod.where((ds_mod.time >= T_start) & (ds_mod.time <= T_end), drop=True)

    # Common time step
    def agg_time(ds_IN, dt_eval, dt_eval_hr, percent_nan_allowed):

        dt_in = (ds_IN.time.values[2] - ds_IN.time.values[1]).astype('timedelta64[h]').astype(float)

        if ((dt_in != dt_eval_hr[dt_eval]) & (
            dt_eval != 'exact')):  # Check if we need to agg (exact skips agg (i.e. snow surveys))
            print('Resampling')

            # Aggregate
            ds_OUT = ds_IN.resample(freq=dt_eval, dim='time', how='mean', label='left', skipna=skipna)

            # Aggregate booleans of not missing, to get fraction in agg period not missing
            obs_fraction_OK = ds_IN.notnull().resample(freq=dt_eval, dim='time', how='mean', label='left')
            if 'p' in ds_IN.data_vars:
                ds_OUT['p'] = ds_IN['p'].resample(freq=dt_eval, dim='time', how='sum', label='left', skipna=skipna)

            # Drop values where fraction exceeds the threshold
            ds_OUT = ds_OUT.where(obs_fraction_OK >= (1-percent_nan_allowed/100) )

        else:
            ds_OUT = ds_IN
        return ds_OUT

    obs_dt_val = agg_time(ds_obs, dt_eval, dt_eval_hr, percent_nan_allowed)
    mod_dt_val = agg_time(ds_mod, dt_eval, dt_eval_hr, percent_nan_allowed)

    # Common non-Missing data (optional)
    if remove_missing:
        print("Removing missing observed timesteps from model")
        I_not_null = (mod_dt_val.notnull()) & (obs_dt_val.notnull())
        mod_dt_val = mod_dt_val.where(I_not_null)
        obs_dt_val = obs_dt_val.where(I_not_null)

    return (obs_dt_val, mod_dt_val)