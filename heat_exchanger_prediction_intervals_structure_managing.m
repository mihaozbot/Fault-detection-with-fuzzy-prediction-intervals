
%If there is only one cluster skip this part
if c>1

    %Dispaly clusters before the evolving mechanisms
    if enable_display_clusters
        heat_exchanger_prediction_intervals_clusters_display
    end

    %Removal mechanism
    if enable_remove
        heat_exchanger_prediction_intervals_removing
    end

    %Merging mechanisms
    merging = 1;
    while merging && (c > 1)

        heat_exchanger_prediction_intervals_merging;

    end

    %Splitting mechanism
    if enable_split
        heat_exchanger_prediction_intervals_splitting
    end

    %Display cluster after the evolving mechanisms
    if enable_display_clusters
        heat_exchanger_prediction_intervals_clusters_display
    end

end

