function Biomechanical_analysis_example()
% An example script to calculate results of ligament tensile
% tests and print them to an excel file.
%
% Aapo Ristaniemi
% Department of Applied Physics, Univeristy of Eastern Finland, 2020

clear all

%----------------------------------------
% Instructions

% 1. Fill in the "User input" below
% 2. Run the script
% 3. Choose the excel-file containing the raw data when asked
% 4. Wait, the results are being calculated
% 5. See the printed figures and Results.xlsx

%----------------------------------------
% User input

% Plot figures?

% type   'no'    to not plot anything
% type   'yes'   to plot all

plot_figures='yes';

% Print results to an excel sheet?

% type   'no'    to not print results
% type   'yes'   to print results to an excel sheet

print_to_excel='yes';

% Use true strain and stress?

% type   'no'    to use engineering strain and stress
% type   'yes'   to use true (logarithmic) strain and stress

true_strain_and_stress='yes';

% End of user input
%----------------------------------------


%----------------------------------------
% Do not change code after this
%----------------------------------------


%----------------------------------------
% Get the raw data file to analyze
[filename,pathname] = uigetfile('*.xlsx');
cd(pathname)
filename((end-4):end)=[];
ligament=filename;

% Some criterions used in the calculations
plot_no=1;
linear_region_criterion=0.08; % 0.08 in the publication
UTS_start_criterion=0.2;     % 0.2 in the publication

% Calculating stress-relaxation results
[sample_length,thickness,width,area,P2E1,P2E2,P2E3,P2E4,zero_mean] = Relaxation_test_analysis(ligament,plot_no,plot_figures,true_strain_and_stress);
plot_no=plot_no+1; % Because 1 figure is printed in relaxation

% Calculating the sinusoidal results
[phase_deg1,phase_deg2,phase_deg3,E_dyn1_MPa,E_dyn2_MPa,E_dyn3_MPa] = Sinusoidal_test_analysis(ligament,sample_length,area,zero_mean,plot_no,plot_figures,true_strain_and_stress);
plot_no=plot_no+1; % Because 1 figure is printed in sinusoidal

% Calculating ultimate tensile test results
[Youngs_modulus_MPa,A_MPa,B_MPa,C_MPa,D_MPa,F,linear_region_length_p,eps_toe_p,sigma_toe_MPa,eps_yield_p,sigma_yield_MPa,Toughness_yield_MPa,strain_at_ultimate_strength_p,Ultimate_strength_MPa,Toughness_failure_MPa] = Ultimate_tensile_test_analysis(ligament,sample_length,area,zero_mean,plot_no,plot_figures,true_strain_and_stress,UTS_start_criterion,linear_region_criterion);

% Printing results to an excel file
if strcmp(print_to_excel,'yes')
    
    % Print file headings
    headings={'Ligament' 'Peak-to-equilibrium ratio for 1st step' 'Peak-to-equilibrium ratio for 2nd step' 'Peak-to-equilibrium ratio for 3rd step' 'Peak-to-equilibrium ratio for 4th step' 'Phase difference at 0.1 Hz' 'Phase difference at 0.5 Hz' 'Phase difference at 1 Hz' 'Dynamic modulus at 0.1 Hz' 'Dynamic modulus at 0.5 Hz' 'Dynamic modulus at 1 Hz' 'Youngs modulus' 'A' 'B' 'C' 'D' 'F' 'Linear region length' 'Toe region strain' 'Toe region stress' 'Yield strain' 'Yield stress' 'Toughness at yield' 'Ultimate strain' 'Ultimate strength' 'Toughness at failure'};
    units={'-' '-' '-' '-' '-' 'degrees' 'degrees' 'degrees' 'MPa' 'MPa' 'MPa' 'MPa' 'MPa' 'MPa' 'MPa' 'MPa' '-' '%' '%' 'MPa' '%' 'MPa' 'mJ/mm^3' '%' 'MPa' 'mJ/mm^3'};
    xlswrite('Results.xlsx',headings,'A1:Z1')
    xlswrite('Results.xlsx',units,'A2:Z2')
    
    % Put all results to same vector
    Results=[P2E1 P2E2 P2E3 P2E4 phase_deg1 phase_deg2 phase_deg3 E_dyn1_MPa E_dyn2_MPa E_dyn3_MPa Youngs_modulus_MPa A_MPa B_MPa C_MPa D_MPa F linear_region_length_p eps_toe_p sigma_toe_MPa eps_yield_p sigma_yield_MPa Toughness_yield_MPa strain_at_ultimate_strength_p Ultimate_strength_MPa Toughness_failure_MPa];
    
    % Identify next free row
    [~,result_file_texts]=xlsread('Results.xlsx');
    printing_row=length(result_file_texts(:,1))+1;

    % Print to excel
    sarake1=['A' num2str(printing_row) ':A' num2str(printing_row)];
    sarake2=['B' num2str(printing_row) ':Z' num2str(printing_row)];
    ligament_r={ligament};
    xlswrite('Results.xlsx',ligament_r,sarake1)
    xlswrite('Results.xlsx',Results,sarake2)

end
end


function [sample_length,thickness,width,area,P2E1,P2E2,P2E3,P2E4,zero_mean] = Relaxation_test_analysis(ligament,plot_no,plot_figures,true_strain_and_stress)
%----------------------------------------
% Stress-relaxation results

excelfilename=[ligament '.xlsx'];
ltw_sheet=[ligament '_l_t_w'];
relax_sheet=[ligament '_relax'];
    
% Get length, thickness and width, and calculate area
l_t_w=xlsread(excelfilename,ltw_sheet);
sample_length=l_t_w(1,1);       % In mm
thickness=l_t_w(1,2)*10^(-3);   % In m
width=l_t_w(1,3)*10^(-3);       % In m
area=(pi/4)*thickness*width;    % NOTE: elliptical assumption. in m^2

% Relaxation steps according to the experiment

step=sample_length*0.02; % in mm
relaxation_steps=[0 step 2*step 3*step 4*step];
 
% Get relaxation data

relax_data=xlsread(excelfilename,relax_sheet);
time_relax=relax_data(:,1);
disp_relax=relax_data(:,2);
grams_relax=relax_data(:,3);

% Now the data is in time-, disp- and grams -vectors

% Subtracting the zero signal (mean of the beginning, where disp is 0)
k=1;
for iiiiii=1:500%length(disp_stressrelax)
    if abs(disp_relax(iiiiii))<0.4
        grams_relax_zero(k)=grams_relax(iiiiii);
        k=k+1;
    end
end

zero_mean=mean(grams_relax_zero);
grams_relax=grams_relax-zero_mean;

% Obtaining time in seconds, disp in mm (tension positive) and force in N

force_relax=9.81*grams_relax/1000;
disp_relax=-disp_relax/1000;    % Disp in mm

%----------------------------------
% Getting starting_row and ending_row of steps
%----------------------------------

nb_steps=length(relaxation_steps);
m=1;

for step_nb=1:(nb_steps-1)

    % Picking up index of step start
    
    for l=m:length(disp_relax)
        difference=abs(disp_relax(l)-relaxation_steps(step_nb));
        if difference>=0.001
            starting_row=l;
            break
        end
    end

    % Picking up index of step end

    for m=starting_row:1:length(disp_relax)
        difference2=abs(disp_relax(m)-relaxation_steps(step_nb+1));
        if difference2<0.001
            ending_row=m;
            break
        end
    end

    % Picking the equilibrium points
    
    equilibrium_point_row(step_nb)=starting_row-1;
    
    % Picking the peak points
    
    peak_point_row(step_nb)=ending_row;
    
end

equilibrium_point_row(end+1)=length(disp_relax);
peak_point_row(end+1)=length(disp_relax); % Just to get equal lengths!

% Getting peak and equilibrium values

for equil_pt_nb=1:length(equilibrium_point_row)
    row_nb=equilibrium_point_row(equil_pt_nb);
    peak_row_nb=peak_point_row(equil_pt_nb);
    epsilon_equilibrium(equil_pt_nb)=disp_relax(row_nb)/sample_length;
    
    sigma_equilibrium(equil_pt_nb)=mean(force_relax((row_nb-50):1:(row_nb)))/area;
    force_equilibrium(equil_pt_nb)=mean(force_relax((row_nb-50):1:(row_nb)));    
    
    force_peak(equil_pt_nb)=force_relax(peak_row_nb);
end

force_equilibrium2345=force_equilibrium;
force_peak1234=force_peak;
sigma_peak1234_MPa=(force_peak/area)/1000000;

% If true stress and strain
if strcmp(true_strain_and_stress,'yes')
    sigma_peak1234_MPa=sigma_peak1234_MPa(1:4).*(1+epsilon_equilibrium(2:end));
    sigma_equilibrium=sigma_equilibrium.*(1+epsilon_equilibrium);
end

%---------------------------------------------------------------
% Calculating peak to equilibrium
sigma_equilibrium_MPa=sigma_equilibrium/1000000;

P2E1=sigma_peak1234_MPa(1)/sigma_equilibrium_MPa(2);
P2E2=sigma_peak1234_MPa(2)/sigma_equilibrium_MPa(3);
P2E3=sigma_peak1234_MPa(3)/sigma_equilibrium_MPa(4);
P2E4=sigma_peak1234_MPa(4)/sigma_equilibrium_MPa(5);

%---------------------------------------------------------------
% These can also be extracted

% f_equil1=force_equilibrium2345(2);
% f_equil2=force_equilibrium2345(3);
% f_equil3=force_equilibrium2345(4);
% f_equil4=force_equilibrium2345(5);
% 
% f_peak1=force_peak1234(1);
% f_peak2=force_peak1234(2);
% f_peak3=force_peak1234(3);
% f_peak4=force_peak1234(4);
% 
% eps_equil1_p=epsilon_equilibrium(2)*100;
% eps_equil2_p=epsilon_equilibrium(3)*100;
% eps_equil3_p=epsilon_equilibrium(4)*100;
% eps_equil4_p=epsilon_equilibrium(5)*100;
% 
% sig_equil1_MPa=sigma_equilibrium_MPa(2);
% sig_equil2_MPa=sigma_equilibrium_MPa(3);
% sig_equil3_MPa=sigma_equilibrium_MPa(4);
% sig_equil4_MPa=sigma_equilibrium_MPa(5);
% 
% sig_peak1_MPa=sigma_peak1234_MPa(1);
% sig_peak2_MPa=sigma_peak1234_MPa(2);
% sig_peak3_MPa=sigma_peak1234_MPa(3);
% sig_peak4_MPa=sigma_peak1234_MPa(4);

stress_relaxation=force_relax/area/1000000;
strain_relaxation=disp_relax/sample_length;

% Calculating true stress if wanted

if strcmp(true_strain_and_stress,'yes')
    stress_relaxation=stress_relaxation.*(1+strain_relaxation);
end

%---------------------------------------------------------------
%  ______ _____ _____ _    _ _____  ______  _____ 
% |  ____|_   _/ ____| |  | |  __ \|  ____|/ ____|
% | |__    | || |  __| |  | | |__) | |__  | (___  
% |  __|   | || | |_ | |  | |  _  /|  __|  \___ \ 
% | |     _| || |__| | |__| | | \ \| |____ ____) |
% |_|    |_____\_____|\____/|_|  \_\______|_____/ 
%---------------------------------------------------------------

% Printing results to a figure

if strcmp(plot_figures,'yes')

    figure(plot_no)

    subplot(3,1,1)
    plot(time_relax,disp_relax)
    xlabel('Time (s)')
    ylabel('Displacement (mm)')
    title([ligament(1:2) ' ' ligament(4:end) ' relaxation test'])
    
    subplot(3,1,2)
    grid on
    plot(time_relax,force_relax)
    xlabel('Time (s)')
    ylabel('Force (N)')

    subplot(3,1,3)
    grid on
    plot(time_relax,stress_relaxation)
    xlabel('Time (s)')
    ylabel('Stress (MPa)')
end
end

function [phase_deg1,phase_deg2,phase_deg3,E_dyn1_MPa,E_dyn2_MPa,E_dyn3_MPa] = Sinusoidal_test_analysis(ligament,sample_length,area,zero_mean,plot_no,plot_figures,true_strain_and_stress)
%----------------------------------------
% Sinusoidal results

excelfilename=[ligament '.xlsx'];
sin1_sheet=[ligament '_sin1'];
sin2_sheet=[ligament '_sin2'];
sin3_sheet=[ligament '_sin3'];

% Getting sinusoidal datas
sin1_data=xlsread(excelfilename,sin1_sheet);
time_sin1=sin1_data(:,1);
disp_sin1=sin1_data(:,2);
grams_sin1=sin1_data(:,3);

% Getting sinusoidal datas
sin2_data=xlsread(excelfilename,sin2_sheet);
time_sin2=sin2_data(:,1);
disp_sin2=sin2_data(:,2);
grams_sin2=sin2_data(:,3);

% Getting sinusoidal datas
sin3_data=xlsread(excelfilename,sin3_sheet);
time_sin3=sin3_data(:,1);
disp_sin3=sin3_data(:,2);
grams_sin3=sin3_data(:,3);

% Now the data is and in time-, disp- and grams -vectors

% Subtracting zero signal

grams_sin1=grams_sin1-zero_mean;
grams_sin2=grams_sin2-zero_mean;
grams_sin3=grams_sin3-zero_mean;

% Obtaining displacement in mm (tension positive) and force in N

force_sin1=9.81*grams_sin1/1000;
disp_sin1=-disp_sin1/1000;    % Disp in mm
force_sin2=9.81*grams_sin2/1000;
disp_sin2=-disp_sin2/1000;    % Disp in mm
force_sin3=9.81*grams_sin3/1000;
disp_sin3=-disp_sin3/1000;    % Disp in mm

% Making correct time vectors
clear ii
time_sin1_2(1)=time_sin1(1);
for ii=2:length(time_sin1)
    time_sin1_2(ii)=time_sin1_2(ii-1)+time_sin1(ii);
end
time_sin1=time_sin1_2'/1000;

clear ii
time_sin2_2(1)=time_sin2(1);
for ii=2:length(time_sin2)
    time_sin2_2(ii)=time_sin2_2(ii-1)+time_sin2(ii);
end
time_sin2=time_sin2_2'/1000;

clear ii
time_sin3_2(1)=time_sin3(1);
for ii=2:length(time_sin3)
    time_sin3_2(ii)=time_sin3_2(ii-1)+time_sin3(ii);
end
time_sin3=time_sin3_2'/1000;

% Calculating stresses and strains

sigma_sin1=force_sin1/area;
epsilon_sin1=disp_sin1/sample_length;
sigma_sin2=force_sin2/area;
epsilon_sin2=disp_sin2/sample_length;
sigma_sin3=(force_sin3/area);
epsilon_sin3=disp_sin3/sample_length;

% Calculating true stress and strain, if wanted

if strcmp(true_strain_and_stress,'yes')
    sigma_sin1=sigma_sin1.*(1+epsilon_sin1);
    epsilon_sin1=log(1+epsilon_sin1);
    
    sigma_sin2=sigma_sin2.*(1+epsilon_sin2);
    epsilon_sin2=log(1+epsilon_sin2);    
    
    sigma_sin3=sigma_sin3.*(1+epsilon_sin3);
    epsilon_sin3=log(1+epsilon_sin3);
end

% Calculating amplitudes

sigma_amp_sin1=(max(sigma_sin1)-min(sigma_sin1))/2;
epsilon_amp_sin1=(max(epsilon_sin1)-min(epsilon_sin1))/2;
sigma_amp_sin2=(max(sigma_sin2)-min(sigma_sin2))/2;
epsilon_amp_sin2=(max(epsilon_sin2)-min(epsilon_sin2))/2;
sigma_amp_sin3=(max(sigma_sin3)-min(sigma_sin3))/2;
epsilon_amp_sin3=(max(epsilon_sin3)-min(epsilon_sin3))/2;

%----------------------------------------
% Fitting sine functions to stress and strain datas

% Frequencies in Hz
freq1=0.1;
freq2=0.5;
freq3=1;

% Guesses for the stresses
%guess=[amplitude,frequency,phase,mean];
guess1 = [sigma_amp_sin1,freq1,pi/8,mean(sigma_sin1)];
guess2 = [sigma_amp_sin2,freq2,pi/8,mean(sigma_sin2)];
guess3 = [sigma_amp_sin3,freq3,pi/2,mean(sigma_sin3)];

% Fits for the stresses
%[p,R,J,CovB,MSE,ErrorModelInfo,x_values_fit,y_values_fit] = Sin_fit(x_values,y_values,guess);
[p1,R,J,CovB,MSE,ErrorModelInfo,time_sin1_fit_sig,sigma_sin1_fit] = Sin_fit(time_sin1,sigma_sin1,guess1);
[p2,R,J,CovB,MSE,ErrorModelInfo,time_sin2_fit_sig,sigma_sin2_fit] = Sin_fit(time_sin2,sigma_sin2,guess2);
[p3,R,J,CovB,MSE,ErrorModelInfo,time_sin3_fit_sig,sigma_sin3_fit] = Sin_fit(time_sin3,sigma_sin3,guess3);

% Guesses for the strains
%guess=[amplitude,frequency,phase,mean];
guess1 = [epsilon_amp_sin1,freq1,pi/8,mean(epsilon_sin1)];
guess2 = [epsilon_amp_sin2,freq2,pi/2,mean(epsilon_sin2)];
guess3 = [epsilon_amp_sin3,freq3,pi/2,mean(epsilon_sin3)];

% Fits for the strains
%[p,R,J,CovB,MSE,ErrorModelInfo,x_values_fit,y_values_fit] = Sin_fit(x_values,y_values,guess);
[p1_eps,R,J,CovB,MSE,ErrorModelInfo,time_sin1_fit_eps,epsilon_sin1_fit] = Sin_fit(time_sin1,epsilon_sin1,guess1);
[p2_eps,R,J,CovB,MSE,ErrorModelInfo,time_sin2_fit_eps,epsilon_sin2_fit] = Sin_fit(time_sin2,epsilon_sin2,guess2);
[p3_eps,R,J,CovB,MSE,ErrorModelInfo,time_sin3_fit_eps,epsilon_sin3_fit] = Sin_fit(time_sin3,epsilon_sin3,guess3);

% Shifting phase to small value, as it may be a multiplication of 2*pi

if p1(3)<0
    while p1(3)<0
        p1(3)=p1(3)+2*pi;
    end
end
if p1(3)>=2*pi
    while p1(3)>=2*pi
        p1(3)=p1(3)-2*pi;
    end
end

if p1_eps(3)<0
    while p1_eps(3)<0
        p1_eps(3)=p1_eps(3)+2*pi;
    end
end
if p1_eps(3)>=2*pi
    while p1_eps(3)>=2*pi
        p1_eps(3)=p1_eps(3)-2*pi;
    end
end

if p2(3)<0
    while p2(3)<0
        p2(3)=p2(3)+2*pi;
    end
end
if p2(3)>=2*pi
    while p2(3)>=2*pi
        p2(3)=p2(3)-2*pi;
    end
end

if p2_eps(3)<0
    while p2_eps(3)<0
        p2_eps(3)=p2_eps(3)+2*pi;
    end
end
if p2_eps(3)>=2*pi
    while p2_eps(3)>=2*pi
        p2_eps(3)=p2_eps(3)-2*pi;
    end
end

if p3(3)<0
    while p3(3)<0
        p3(3)=p3(3)+2*pi;
    end
end
if p3(3)>=2*pi
    while p3(3)>=2*pi
        p3(3)=p3(3)-2*pi;
    end
end

if p3_eps(3)<0
    while p3_eps(3)<0
        p3_eps(3)=p3_eps(3)+2*pi;
    end
end
if p3_eps(3)>=2*pi
    while p3_eps(3)>=2*pi
        p3_eps(3)=p3_eps(3)-2*pi;
    end
end

% Calculating phase difference
phase_rad1=-p1_eps(3)+p1(3);
phase_deg1=phase_rad1*(360/(2*pi));

phase_rad2=-p2_eps(3)+p2(3);
phase_deg2=phase_rad2*(360/(2*pi));

phase_rad3=-p3_eps(3)+p3(3);
phase_deg3=phase_rad3*(360/(2*pi));

phase_degrees=[phase_deg1,phase_deg2,phase_deg3];

%----------------------------------------
% Calculating dynamic stiffnesses

dynamic_modulus_sin1=p1(1)/p1_eps(1);
dynamic_modulus_sin2=p2(1)/p2_eps(1);
dynamic_modulus_sin3=p3(1)/p3_eps(1);

Dynamic_moduli_MPa=[dynamic_modulus_sin1,dynamic_modulus_sin2,dynamic_modulus_sin3]/1000000;

%----------------------------------------
% Storage modulus and loss modulus can be calculated

% Storage_modulus_sin1=p1(1)/p1_eps(1)*cos(phase_rad1);
% Loss_modulus_sin1=p1(1)/p1_eps(1)*sin(phase_rad1);
% Storage_modulus_sin2=p2(1)/p2_eps(1)*cos(phase_rad2);
% Loss_modulus_sin2=p2(1)/p2_eps(1)*sin(phase_rad2);
% Storage_modulus_sin3=p3(1)/p3_eps(1)*cos(phase_rad3);
% Loss_modulus_sin3=p3(1)/p3_eps(1)*sin(phase_rad3);
% 
% Storage_moduli_MPa=[Storage_modulus_sin1,Storage_modulus_sin2,Storage_modulus_sin3]/1000000;
% Loss_moduli_MPa=[Loss_modulus_sin1,Loss_modulus_sin2,Loss_modulus_sin3]/1000000;

%---------------------------------------------------------------
% Adjusting prints

E_dyn1_MPa=Dynamic_moduli_MPa(1);
E_dyn2_MPa=Dynamic_moduli_MPa(2);
E_dyn3_MPa=Dynamic_moduli_MPa(3);

% E_storage1_MPa=Storage_moduli_MPa(1);
% E_storage2_MPa=Storage_moduli_MPa(2);
% E_storage3_MPa=Storage_moduli_MPa(3);
% 
% E_loss1_MPa=Loss_moduli_MPa(1);
% E_loss2_MPa=Loss_moduli_MPa(2);
% E_loss3_MPa=Loss_moduli_MPa(3);

%---------------------------------------------------------------
%  ______ _____ _____ _    _ _____  ______  _____ 
% |  ____|_   _/ ____| |  | |  __ \|  ____|/ ____|
% | |__    | || |  __| |  | | |__) | |__  | (___  
% |  __|   | || | |_ | |  | |  _  /|  __|  \___ \ 
% | |     _| || |__| | |__| | | \ \| |____ ____) |
% |_|    |_____\_____|\____/|_|  \_\______|_____/ 
%---------------------------------------------------------------

% PRINTING results to figures

% Plotting signals, just to see that they are reasonable

if strcmp(plot_figures,'yes')
    
    figure(plot_no)
    subplot(3,2,1)
    yyaxis left          % plot against left y-axis Requires Matlab R2016!!!
    plot(time_sin1,disp_sin1,'b');%,time_sin1,force_sin1,'r')
    ylabel('Displacement (mm)')
    yyaxis right         % plot against right y-axis
    plot(time_sin1,force_sin1,'r')
    xlabel('Time (s)')
    ylabel('Force (N)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin1, visual check that forces and displacements are OK'])
    
    subplot(3,2,3)
    yyaxis left          % plot against left y-axis  
    plot(time_sin2,disp_sin2,'b')
    ylabel('Displacement (mm)')
    yyaxis right         % plot against right y-axis
    plot(time_sin2,force_sin2,'r')
    xlabel('Time (s)')
    ylabel('Force (N)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin2'])

    subplot(3,2,5)
    yyaxis left          % plot against left y-axis  
    plot(time_sin3,disp_sin3,'b')
    ylabel('Displacement (mm)')
    yyaxis right         % plot against right y-axis
    plot(time_sin3,force_sin3,'r')
    xlabel('Time (s)')
    ylabel('Force (N)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin3'])

    % Stresses and strains, original and fit
    subplot(3,2,2)
    yyaxis left          % plot against left y-axis Requires Matlab R2016!!!
    plot(time_sin1,epsilon_sin1,time_sin1_fit_eps,epsilon_sin1_fit,'b');
    ylabel('Epsilon (-)')
    yyaxis right         % plot against right y-axis
    plot(time_sin1,sigma_sin1,time_sin1_fit_sig,sigma_sin1_fit,'r')
    xlabel('Time (s)')
    ylabel('Stress (Pa)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin1, visual check that fits are OK'])

    subplot(3,2,4)
    yyaxis left          % plot against left y-axis  
    plot(time_sin2,epsilon_sin2,time_sin2_fit_eps,epsilon_sin2_fit,'b')
    ylabel('Epsilon (-)')
    yyaxis right         % plot against right y-axis
    plot(time_sin2,sigma_sin2,time_sin2_fit_sig,sigma_sin2_fit,'r')
    xlabel('Time (s)')
    ylabel('Stress (Pa)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin2'])

    subplot(3,2,6)
    yyaxis left          % plot against left y-axis  
    plot(time_sin3,epsilon_sin3,time_sin3_fit_eps,epsilon_sin3_fit,'b')
    ylabel('Epsilon (-)')
    yyaxis right         % plot against right y-axis
    plot(time_sin3,sigma_sin3,time_sin3_fit_sig,sigma_sin3_fit,'r')
    xlabel('Time (s)')
    ylabel('Stress (Pa)')
    title([ligament(1:2) ' ' ligament(4:end) ' sin3'])    

end
end

function [p,J,R,CovB,MSE,ErrorModelInfo,x_values_fit,y_values_fit] = Sin_fit(x_values,y_values,guess)
%-------------------------------
% Fitting a sinusoidal function to data

% Using nlinfit

% Giving function form
fun = @(p,x_values) p(1)*sin(2*pi*p(2)*x_values+p(3))+p(4);

% Setting options
options=statset('MaxIter',1000);

% Solving everything with nlinfit
[p,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x_values,y_values,fun,guess,options);

% Negative amplitude made positive using properties of sine function :)

if p(1)<0
    p(3)=p(3)+pi;
    p(1)=-p(1);
end
%-------------------------------
% Optional: lsqcurvefit, gives the same result
%[p,fminres] = lsqcurvefit(fun,guess,x_values,y_values);

%-------------------------------
% Curve fit

x_values_fit=x_values(1):0.001:x_values(end);
y_values_fit=(p(1)*sin(2*pi*p(2)*x_values_fit+p(3))+p(4));

end

function[Youngs_modulus_MPa,A_MPa,B_MPa,C_MPa,D_MPa,F,linear_region_length_p,eps_toe_p,sigma_toe_MPa,eps_yield_p,sigma_yield_MPa,Toughness_yield_MPa,strain_at_ultimate_strength_p,Ultimate_strength_MPa,Toughness_failure_MPa] = Ultimate_tensile_test_analysis(ligament,sample_length,area,zero_mean,plot_no,plot_figures,true_strain_and_stress,UTS_start_criterion,linear_region_criterion)
%----------------------------------------------------------
% Function to analyze ultimate tensile test data

% Getting the ultimate data
excelfilename=[ligament '.xlsx'];
UTS_sheet=[ligament '_UTS'];
UTS_data=xlsread(excelfilename,UTS_sheet);
time_UTS=UTS_data(:,1);
disp_UTS=UTS_data(:,2);
grams_UTS=UTS_data(:,3);

% Now the data is in time-, disp- and grams -vectors

grams_UTS=grams_UTS-zero_mean;

% Obtaining time in seconds, disp in mm (tension positive) and force in N

force_UTS=9.81*grams_UTS/1000;
disp_UTS=-disp_UTS/1000;

% Filtering the force signal if needed
sample_freq=100;
[force_UTS] = Filter(sample_freq,time_UTS,force_UTS);

% Ignoring the first 3600 seconds of recovery time
for row_nb=1:length(disp_UTS)
    if abs(time_UTS(row_nb))>3595
        starting_row1=row_nb;
        break
    end
end

for row_nb2=starting_row1:length(disp_UTS)
    if abs(disp_UTS(row_nb2))>0.001
        starting_row_UTS=row_nb2-1; % Real starting row
        break
    end
end

disp_UTS2=disp_UTS(starting_row_UTS:end);
force_UTS2=force_UTS(starting_row_UTS:end);

% Now we have only the data of the test (no 3600 s recovery)

%----------------------------------------------------------
% Original engineering stress and strain
eps_eng=disp_UTS2/sample_length;
sigma_eng=force_UTS2/area;

eps_UTS=eps_eng;
sigma_UTS=sigma_eng;

sample_length2=sample_length;
eps_shift_eng=0;

%----------------------------------------------------------
% Changing to true stress and strain if wanted

if strcmp(true_strain_and_stress,'yes')
    eps_UTS=log(1+eps_eng);
    sigma_UTS=sigma_eng.*(1+eps_eng);
end

%----------------------------------------------------------
% Shifting the data with shift criterion (UTS_start_criterion)

UTS_start_criterion=UTS_start_criterion-0.05;
starting_row_of_UTS_mod=0;
if not(UTS_start_criterion==0)
    for sigma_row=6:length(sigma_UTS) % Going through all the sigmas
    sigma_value=mean(sigma_UTS((sigma_row-5):1:(sigma_row+5))); % Picking the stress value of measured data
        if sigma_value>UTS_start_criterion*1000000 % In Pa
            starting_row_of_UTS_mod=sigma_row;
            break
        end 
    end
    
    sigma_eng2=sigma_eng(starting_row_of_UTS_mod:end);
    sigma_UTS=sigma_eng2;
    
    % Adjusting zero-load length
    disp_UTS_mod=disp_UTS2(starting_row_of_UTS_mod:end);
    sample_length_mod=sample_length+disp_UTS_mod(1);
    
    % Shifting disp to 0
    shift_disp=disp_UTS_mod(1);
    eps_shift_eng=100*shift_disp/sample_length;

    disp_UTS_mod=disp_UTS_mod-disp_UTS_mod(1);
    % Calculating new engineering (and true) strain
    eps_eng2=disp_UTS_mod/sample_length_mod;
    eps_UTS=eps_eng2;
    
    sample_length2=sample_length_mod;
    
    if strcmp(true_strain_and_stress,'yes')
        eps_UTS=log(1+eps_eng2);
        sigma_UTS=sigma_eng2.*(1+eps_eng2);
    end
end

% Now the data is shifted (if start criterion for UTS test is used)
% and is either engineering or true.

%----------------------------------------------------------
% Calculating Young's modulus of the linear region

% Finding the length of interval in terms of row numbers

for index_nb=1:length(eps_UTS)
    if abs(eps_UTS(index_nb))>linear_region_criterion
        nb_rows_in_interval_UTS=index_nb;
        break
    end
end

% Now we have the number of rows in interval

% Making fits and finding the maximum of them

fit_no=1;

for s_row=1:(length(eps_UTS)-nb_rows_in_interval_UTS)
    coeff_UTS=polyfit(eps_UTS(s_row:(s_row+nb_rows_in_interval_UTS)),sigma_UTS(s_row:(s_row+nb_rows_in_interval_UTS)),1);
    moduli(fit_no)=coeff_UTS(1);
    s_row_no(fit_no)=s_row;
    fit_no=fit_no+1;
end

[max_modulus,max_modulus_location]=max(moduli);
starting_row_number=s_row_no(max_modulus_location);

% Obtaining the max modulus fit
coeff_UTS_max=polyfit(eps_UTS(starting_row_number:(starting_row_number+nb_rows_in_interval_UTS)),sigma_UTS(starting_row_number:(starting_row_number+nb_rows_in_interval_UTS)),1);
E_UTS_MPa=coeff_UTS_max(1)/1000000;
Youngs_modulus_MPa=E_UTS_MPa;

%-----------------------------------------------------
% Calculating 2-order fit for the toe region

% Finding the end of toe region (toe region assumed to start at eps=0)
% Criterion: when difference with the max modulus fit is delta_eps=0.006

criterion=0.006;

for sigma_row=1:length(sigma_UTS) % Going through all the sigmas
    sigma_value=sigma_UTS(sigma_row); % Picking the stress value of measured data
    eps_fit=(sigma_value-coeff_UTS_max(2))/coeff_UTS_max(1); % Epsilon value of the fit
    eps_value=eps_UTS(sigma_row); % epsilon of the measured data
    if abs(eps_fit-eps_value)<criterion
        ending_row_of_toe_region=sigma_row;
        break
    end 
end

sigma_toe_region=sigma_UTS(1:ending_row_of_toe_region);
eps_toe_region=eps_UTS(1:ending_row_of_toe_region);

sigma_toe_region_end_MPa=sigma_toe_region(end)/1000000;
eps_toe_region_end_p=eps_toe_region(end)*100;

sigma_toe_MPa=sigma_toe_region_end_MPa;
eps_toe_p=eps_toe_region_end_p;

% Making fit, 2nd order

toe_fun= @(toe_fit,eps_toe_region) toe_fit(1)+toe_fit(2)*eps_toe_region+toe_fit(3)*eps_toe_region.^2;
guess=[1 1 1];

[toe_fit,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(eps_toe_region,sigma_toe_region,toe_fun,guess);

% fit for plotting
eps_toe_fit=0:0.000001:(eps_toe_region(end)+0.02);
sigma_toe_fit=toe_fit(1)+toe_fit(2)*eps_toe_fit+toe_fit(3)*eps_toe_fit.^2;

A_MPa=toe_fit(3)/1000000;
B_MPa=toe_fit(2)/1000000;
C_MPa=toe_fit(1)/1000000;

%-----------------------------------------------
% Making toe region fit according to Fung (1967)

fung_fun= @(fung_fit,eps_toe_region) fung_fit(1)*(exp(fung_fit(2)*eps_toe_region)-1);
fung_guess=[3*10^5 30];
lb=[];
ub=[];
fung_options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',1000);
fung_fit = lsqcurvefit(fung_fun,fung_guess,eps_toe_region,sigma_toe_region,lb,ub,fung_options);

% fit for plotting
eps_fung_fit=0:0.000001:(eps_toe_region(end)+0.02);
sigma_fung_fit=fung_fit(1)*(exp(fung_fit(2)*eps_fung_fit)-1);

D_MPa=fung_fit(1)/1000000;
F=fung_fit(2);

%---------------------------------
% Calculating yield point

% Finding the yield point
% Criterion: when difference with the max modulus fit is 0.006.

yield_criterion=0.006;

for sigma_row2=(ending_row_of_toe_region+round(nb_rows_in_interval_UTS/2)):length(sigma_UTS) % Going through all the sigmas after toe region
    sigma_value=sigma_UTS(sigma_row2); % Picking the stress value of measured data
    eps_fit=(sigma_value-coeff_UTS_max(2))/coeff_UTS_max(1); % Epsilon value of the fit
    eps_value=eps_UTS(sigma_row2); % epsilon of the measured data
    if abs(eps_fit-eps_value)>yield_criterion
        yield_row=sigma_row2;
        break
    end 
end

% Yield point
sigma_yield_MPa=sigma_UTS(yield_row)/1000000;
eps_yield_p=eps_UTS(yield_row)*100;

% End of yield point calculation
%---------------------------------

% Calculating fit for the linear region
max_modulus_eps=[(eps_toe_region(end)-0.013):0.001:(eps_yield_p/100+0.013)];
max_modulus_sig_MPa=(coeff_UTS_max(1)*max_modulus_eps+coeff_UTS_max(2))/1000000;

max_modulus_eps_p=max_modulus_eps*100; % For printing

% Linear region length
linear_region_length_p=(eps_UTS(yield_row)-eps_UTS(ending_row_of_toe_region))*100;

%-------------------------------------------
% Ultimate strength

[Ultimate_strength,ultimate_location]=max(sigma_UTS);
[Ultimate_force,ultimate_location_force]=max(force_UTS2);
Ultimate_strength_MPa=Ultimate_strength/1000000;

% Strain at ultimate strength

strain_at_ultimate_strength_p=eps_UTS(ultimate_location)*100;

%-------------------------------------------
% Toughness at failure

eps_uts_area=eps_UTS(1:ultimate_location);
sigma_uts_area=sigma_UTS(1:ultimate_location);
sigma_uts_area(end+1)=0;
sigma_uts_area(end+1)=0;
eps_uts_area(end+1)=eps_uts_area(end);
eps_uts_area(end+1)=eps_uts_area(1);
Toughness_failure_MPa=polyarea(eps_uts_area,sigma_uts_area)/1000000;

%-------------------------------------------
% Toughness at yield

eps_yield_area=eps_UTS(1:yield_row);
sigma_yield_area=sigma_UTS(1:yield_row);
sigma_yield_area(end+1)=0;
sigma_yield_area(end+1)=0;
eps_yield_area(end+1)=eps_yield_area(end);
eps_yield_area(end+1)=eps_yield_area(1);
Toughness_yield_MPa=polyarea(eps_yield_area,sigma_yield_area)/1000000;

%---------------------------------------------------------------
%  ______ _____ _____ _    _ _____  ______  _____ 
% |  ____|_   _/ ____| |  | |  __ \|  ____|/ ____|
% | |__    | || |  __| |  | | |__) | |__  | (___  
% |  __|   | || | |_ | |  | |  _  /|  __|  \___ \ 
% | |     _| || |__| | |__| | | \ \| |____ ____) |
% |_|    |_____\_____|\____/|_|  \_\______|_____/ 
%---------------------------------------------------------------

if strcmp(plot_figures,'yes')

    figure(plot_no)
    grid on

    plot(eps_UTS*100,sigma_UTS/1000000,max_modulus_eps_p,max_modulus_sig_MPa,eps_toe_region*100,sigma_toe_region/1000000,'--',eps_toe_region_end_p,sigma_toe_region_end_MPa,'o',eps_toe_fit*100,sigma_toe_fit/1000000,eps_fung_fit*100,sigma_fung_fit/1000000,eps_yield_p,sigma_yield_MPa,'o',strain_at_ultimate_strength_p,Ultimate_strength_MPa,'o')
    xlabel('Strain (%)')
    ylabel('Stress (MPa)')
    legend(['UTS test'],['E_l_i_n_e_a_r= ' num2str(E_UTS_MPa) ' MPa'],'Toe region',['Toe region end, eps=' num2str(eps_toe_region_end_p) ' %, sig=' num2str(sigma_toe_region_end_MPa) ' MPa'],'Toe fit','Fung fit',['Yield point, eps=' num2str(eps_yield_p) ' %, sig=' num2str(sigma_yield_MPa) ' MPa'],['Ultimate point, eps=' num2str(strain_at_ultimate_strength_p) ' %, sig=' num2str(Ultimate_strength_MPa) ' MPa'],'location','southeast')
    title([ligament(1:2) ' ' ligament(4:end) ' UTS test'])

end
end

function [low_data] = Filter(sample_freq,time,force)
% Filtering out noise from the signal

% Butterworth filter

data=force; 
f=sample_freq;%  sampling frequency
f_cutoff = 10; % cutoff frequency

fnorm = f_cutoff/(f/2); % normalized cut off freq, you can change it to any value depending on your requirements

[b1,a1] = butter(4,fnorm,'low'); % Low pass Butterworth filter of order 4
low_data = filtfilt(b1,a1,data); % filtering

end
