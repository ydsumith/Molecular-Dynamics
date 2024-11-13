% Conversion factors
eV_to_kcal_per_mol = 23.0605;
g_per_cc_to_kg_per_m3 = 1000;

% Initialize variables
file_indices = 250:10:650;
num_files = length(file_indices);

% Excel filenames
density_excel_filename = 'average_density_profiles_kg_per_m3.xlsx';
pe_excel_filename = 'average_pe_profiles_kcal_per_mol.xlsx';
ke_excel_filename = 'average_ke_profiles_kcal_per_mol.xlsx';
pressure_excel_filename = 'average_pressure_profiles_bar.xlsx';

% Delete existing Excel files if they exist
if exist(density_excel_filename, 'file') == 2
    delete(density_excel_filename);
end
if exist(pe_excel_filename, 'file') == 2
    delete(pe_excel_filename);
end
if exist(ke_excel_filename, 'file') == 2
    delete(ke_excel_filename);
end
if exist(pressure_excel_filename, 'file') == 2
    delete(pressure_excel_filename);
end

% Process each profile type and write to separate files
write_profiles(file_indices, 'density_profile_', 'Density (kg/m^3)', g_per_cc_to_kg_per_m3, density_excel_filename);
write_profiles(file_indices, 'pe_profile_', 'PE (kcal/mol)', eV_to_kcal_per_mol, pe_excel_filename);
write_profiles(file_indices, 'ke_profile_', 'KE (kcal/mol)', eV_to_kcal_per_mol, ke_excel_filename);
write_pressure_profiles(file_indices, 'pressure_profile_', 'Pressure (bar)', pressure_excel_filename);

% Display a message confirming the export
fprintf('Data from all files have been written to separate Excel files:\n');
fprintf(' - %s\n', density_excel_filename);
fprintf(' - %s\n', pe_excel_filename);
fprintf(' - %s\n', ke_excel_filename);
fprintf(' - %s\n', pressure_excel_filename);

% Function to process profiles and write to an Excel file
function write_profiles(file_indices, file_prefix, label, conversion_factor, excel_filename)
    num_files = length(file_indices);

    % Initialize Excel data container
    excel_data = {};

    % Add header row
    header = {'Coord1'};
    for idx = 1:num_files
        file_index = file_indices(idx);
        header{end+1} = sprintf('%s_%d', label, file_index);
    end
    excel_data = [excel_data; header];

    % Initialize container for Coord1 and profile values
    coord1_data = [];
    profile_values = [];

    for idx = 1:num_files
        file_index = file_indices(idx);
        filename = sprintf('%s%d.txt', file_prefix, file_index);
        data_table = read_profile_file(filename, conversion_factor);

        % Store Coord1 once
        if isempty(coord1_data)
            coord1_data = data_table.Coord1;
        end

        % Add the computed averages to profile values
        profile_values(:, idx) = data_table.Average_Value; %#ok<AGROW>
    end

    % Combine Coord1 and profile values
    excel_data = [excel_data; num2cell([coord1_data, profile_values])];

    % Write data to Excel
    writecell(excel_data, excel_filename);
end

% Function to process pressure profiles and write to an Excel file
function write_pressure_profiles(file_indices, file_prefix, label, excel_filename)
    num_files = length(file_indices);

    % Initialize Excel data container
    excel_data = {};

    % Add header row
    header = {'Coord1'};
    for idx = 1:num_files
        file_index = file_indices(idx);
        header{end+1} = sprintf('%s_%d', label, file_index);
    end
    excel_data = [excel_data; header];

    % Initialize container for Coord1 and profile values
    coord1_data = [];
    pressure_values = [];

    for idx = 1:num_files
        file_index = file_indices(idx);
        filename = sprintf('%s%d.txt', file_prefix, file_index);
        data_table = read_pressure_profile(filename);

        % Store Coord1 once
        if isempty(coord1_data)
            coord1_data = data_table.Coord1;
        end

        % Add the computed averages to pressure values
        pressure_values(:, idx) = data_table.Average_Pressure; %#ok<AGROW>
    end

    % Combine Coord1 and pressure values
    excel_data = [excel_data; num2cell([coord1_data, pressure_values])];

    % Write data to Excel
    writecell(excel_data, excel_filename);
end

% Function to read profile files and calculate averages with conversion
function data_table = read_profile_file(filename, conversion_factor)
    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Skip header lines
    fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);

    % Read the data
    data = textscan(fid, '%f %f %f %f');
    fclose(fid);

    % Extract relevant columns
    coord1 = round(data{2}, 5); % Coord1 rounded to avoid precision issues
    value = data{4} * conversion_factor; % Apply conversion factor

    % Group data by Coord1 and compute averages
    [unique_coord1, ~, idx_coord1] = unique(coord1);
    avg_value = accumarray(idx_coord1, value, [], @mean);

    % Create a table
    data_table = table(unique_coord1, avg_value, 'VariableNames', {'Coord1', 'Average_Value'});
end

% Function to read pressure profile files and calculate averages
function data_table = read_pressure_profile(filename)
    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Skip header lines
    fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);

    % Read the data
    data = textscan(fid, '%f %f %f %f %f %f');
    fclose(fid);

    % Extract relevant columns
    coord1 = round(data{2}, 5); % Coord1 rounded to avoid precision issues
    sxx = data{4};
    syy = data{5};
    szz = data{6};

    % Compute pressure
    %P_conversion_factor = 1602.1766 * 1e4; % eV/Å³ to bars
    %P_conversion_factor = 1 / 1.873156418e+7; % normalizing with volume total 
    P_conversion_factor = 1 / 24183.83192974265; % normalizing with volume slab
    pressure = -P_conversion_factor*(sxx + syy + szz) / 3;

    % Filter out zero values
    non_zero_idx = pressure ~= 0;
    coord1 = coord1(non_zero_idx);
    pressure = pressure(non_zero_idx);

    % Group data by Coord1 and compute averages
    [unique_coord1, ~, idx_coord1] = unique(coord1);
    avg_pressure = accumarray(idx_coord1, pressure, [], @mean);

    % Create a table
    data_table = table(unique_coord1, avg_pressure, 'VariableNames', {'Coord1', 'Average_Pressure'});
end
