def newvalue_fun(time, model, netcdf_path_format, is_historical=False, scenario=None):
    newvalue_dict = {}
    projections_dict = {}

    for model_name in model_names:
        for species_name in species_names:
            if is_historical:
                ds = xr.open_dataset(netcdf_path_format.format(model, taxa, species_name, time), decode_times=False)
            else:
                ds = xr.open_dataset(netcdf_path_format.format(model, taxa, model_name, scenario, species_name, time), decode_times=False)

            newvalue = ds["newvalue"]
            newvalue_dict[(model, model_name, species_name)] = newvalue

    for model_name in model_names:
        value_list = []
        for species_name in species_names:
            value_bin = newvalue_dict[(model, model_name, species_name)]
            value_list.append(value_bin)


        value_bin_list = list(value_list)
        mean_value_bin = xr.concat(value_bin_list, dim="species").sum(dim="species")  # Sum over species
        projections_dict[model_name] = mean_value_bin

    return mean_value_bin
