% Miha Ožbot
% Final version: 28.11.2022
% Reviewed: 06.11.2023
clc; clear all; close all;
set(0,'defaultTextInterpreter','latex')
warning('off','all')

%Enable flags, if set to true (or any number) activates the computation of
% the related code segment
enable_identification = 1; %Enable the identification of model, else load
enable_display_clusters = 0; %Display clusters during identification
enable_split = 0; %Enable splitting mechanism
enable_remove = 1; %Enable removal mechanism
enable_disturbance = 1; %Enable the use of the reservoir signal in the regressor
enable_synthetic_fault = 0; % Enable additional faults 

if enable_identification

    %Parameters
    u_min = 4; %Min input signal
    u_max = 20; %Max input signal
    k_split = 0.85; % Principal component split value
    kappa_K = 0.3; %Merge dissimilarity measure value
    kappa_olap = 6; %Cluster overlapping merging threshold
    kappa_r = 0.1; %affine constant merging threshold
    kappa_split = 0; %Rule splitting threshold, if 0 a rule is always split
    kappa_remove = 100^2; %Rule removal threshold

    kappa_c = 0.01; %Minimal confidence interval for adaptation
    m = 1; %Number of delayed inputs used in the model
    n = 1; %Number of delayed outputs used in the model
    P_alpha = 1e10; %heuristic parameter for the information matrix
    etta = 0.5; %Fuzziness factor
    A_f = zeros(1,n+1); A_f(1) = 1; %Filter for the optimization
    t_deviations = 2; %Number of standard deviations for the confidence interval
    Gamma_min = exp(-2^2); %Minimal activation required to adapt a rule
    S_0 = 1e1; %Minimal variance of clusters

    %Initialize the model
    heat_exchanger_prediction_intervals_initialization

    %Single-pass over the data
    for k = k_0:1:length(u)

        %Data acquisition for measurements or data
        heat_exchanger_prediction_intervals_data_acquisition

        %Compute distance
        Gamma = zeros(c,1);
        d2 = zeros(c,1);
        Sigma = zeros(n_z,n_z,c);
        delta = zeros(2, 2, c);

        if c > 0
            for i = 1:c

                %Calculation of clusters projection parameter
                Sigma(:,:,i) = S_c(:,:,i)/n_c(i);
                [~,S,V] = svd(Sigma(:,:,i));
                delta(:,:,i) = diag([max(abs(V(1,:)'.*diag(sqrt(S)))),...
                    max(abs(V(2,:)'.*diag(sqrt(S))))]);

                %d2(i) = (z(k,1:2) - mi(i,1:2))*pinv(delta(:,:,i))*(z(k,1:2)- mi(i,1:2))';
                d2(i) = (z(k,1) - mi(i,1))*pinv(delta(1,1,i))*(z(k,1)- mi(i,1))';

                if d2(i)<0
                    disp('Error distance!')
                    d2(i) = Inf;
                end

                Gamma(i) = exp(-d2(i)^etta) + 1e-15;

            end
        end

        %Computer maximum activation
        [~,jg] = max(Gamma);
        NGamma = Gamma/sum(Gamma);

        if (Gamma(jg) > Gamma_min) + ~(abs(u(k-1)-u(k)) > 0)

            %Incremental identification of the parameters and structure
            heat_exchanger_prediction_intervals_optimization

            %Adapt the filter used in the parameter optimization
            heat_exchanger_prediction_intervals_filter_adaptation

            %Compute the mean square error of the model
            heat_exchanger_prediction_intervals_error_accumulation

        else

            %The evolving mechanisms, merging, splitting, removoval of rules
            heat_exchanger_prediction_intervals_structure_managing

            %Add new rule
            heat_exchanger_prediction_intervals_rule_addition

        end

        %Number of clusters at time step k for review
        c_all(k) = c;

        %Dispaly clusters when input changes
        if abs(u(k-1)-u(k)) > 0

            heat_exchanger_prediction_intervals_verification
            heat_exchanger_prediction_intervals_clusters_display
        end

    end

    %The evolving mechanisms, merging, splitting, removoval of rules
    heat_exchanger_prediction_intervals_structure_managing

    %Save the model and code state
    save heat_exchanger_prediction_intervals_identification

end

%load model
load heat_exchanger_prediction_intervals_identification

%Dispaly clusters
heat_exchanger_prediction_intervals_clusters_display

%Verification of the identified model
heat_exchanger_prediction_intervals_verification

%Fault detection with N-steps ahead prediction and prediction intervals
fault_threshold = 0.002; %Threshold for alarm
k_ahead = 10; %Number prediction steps ahead
heat_exchanger_prediction_intervals_n_step_prediction

%Fault detection with model simulation and confidence itervals
%heat_exchanger_prediction_intervals_simulation
