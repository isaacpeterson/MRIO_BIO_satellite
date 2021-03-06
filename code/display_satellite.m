
function display_satellite(satellite_object, satellite_params, species_characteristics, sector_lengths, sorted_country_names)

    ytick_object = build_species_labels(species_characteristics.sorted_species_classification, satellite_params.species_to_label);
    
    if strcmp(satellite_params.display_type, 'global')
        xtick_object = build_country_labels(satellite_params.sectors_to_label, sector_lengths, sorted_country_names);
    else 
        xtick_object = build_sector_labels(satellite_params.sectors_to_label, sector_lengths, sorted_country_names);
    end
    
    if (satellite_params.display_domestic_satellite == true)
        figure(1)
        display_current_satellite([satellite_object.domestic_satellite{:}], 'domestic threats satellite', ytick_object, xtick_object)
    end
    
    if (satellite_params.display_global_satellite == true)
        figure(2)
        display_current_satellite([satellite_object.global_satellite{:}], 'global threats satellite', ytick_object, xtick_object)
    end
    
    if (satellite_params.display_global_satellite == true)
        figure(3)
        display_current_satellite([satellite_object.domestic_satellite{:}] + [satellite_object.global_satellite{:}], 'total threats satellite', ytick_object, xtick_object)
    end
    
end


function display_current_satellite(current_satellite, plot_title, ytick_object, xtick_object)
    colormap hot
    imagesc(log(current_satellite)) 
    set(gca, 'ytick', ytick_object.species_ticks_to_use)
    set(gca,'YtickLabel', ytick_object.species_labels_to_use)
    set(gca, 'xtick', xtick_object.sector_ticks)
    set(gca, 'XtickLabel', xtick_object.sector_labels)
    set(gca, 'XTickLabelRotation', 90)
    title(plot_title)
end

function tick_object = build_sector_labels(countries_to_label, sector_lengths, sorted_sector_names)
    sector_vec = build_sector_vec(sector_lengths);
    if strcmp(countries_to_label, 'all')
        sectors_to_label = find(sector_lengths > 0);
    else 
        sectors_to_label = find(ismember(sorted_sector_names, countries_to_label));
    end   
    
    tick_object = struct();
    tick_object.sector_labels = sorted_sector_names(sectors_to_label);
    sector_vec_to_use = [0, sector_vec];
    tick_object.sector_ticks = floor(0.5*(sector_vec_to_use(sectors_to_label) + sector_vec_to_use(sectors_to_label + 1)));
end

function tick_object = build_country_labels(countries_to_label, sector_lengths, sorted_country_names)
    sector_vec = build_sector_vec(sector_lengths);
    if strcmp(countries_to_label, 'all')
        sectors_to_label = find(sector_lengths > 0);
    else 
        sectors_to_label = find(ismember(sorted_country_names, countries_to_label));
    end   
    
    tick_object = struct();
    tick_object.sector_labels = sorted_country_names(sectors_to_label);
    sector_vec_to_use = [0, sector_vec];
    tick_object.sector_ticks = floor(0.5*(sector_vec_to_use(sectors_to_label) + sector_vec_to_use(sectors_to_label + 1)));
end

function tick_object = build_species_labels(classes_to_use, species_to_label)
    
    [species_labels, species_class_counter] = unique(classes_to_use, 'stable');
    if strcmp(species_to_label, 'all')
        species_to_label = species_labels;
    end
    [~, ~, species_inds_to_use] = intersect(species_to_label, species_labels);
    species_class_counter = vertcat(species_class_counter, numel(classes_to_use));
    tick_object = struct();
    tick_object.species_ticks_to_use = floor(0.5*(species_class_counter(species_inds_to_use) + species_class_counter(species_inds_to_use + 1)));
    tick_object.species_labels_to_use = species_labels(species_inds_to_use);
    
end


function [sector_vec] = build_sector_vec(sector_lengths)
    current_NCOUN = length(sector_lengths);
    sector_vec = zeros(1, current_NCOUN);
    for country_ind = 1:current_NCOUN
        sector_vec(country_ind) = sum(sector_lengths(1:(country_ind)));
    end
end