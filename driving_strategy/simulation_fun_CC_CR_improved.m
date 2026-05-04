
% function [Total_E] = simulation_fun(X)
% function [obj] = simulation_fun(X)
% function [running_inter] = simulation_fun(X)
%function [running_inter,Total_E]= simulation_fun(X)
function [running_inter, Total_E, s_out, v_out, vlim_out] = simulation_fun_CC_CR_improved(X)
% Improved driving strategy simulator for CC-CR (Cruise-Coasting-Cruise) mode
% Enhanced version with adaptive acceleration/braking scaling parameters
%
% Inputs:
%   X = decision variables with optional adaptive END parameters:
%       [...base_decision_vars..., alpha_trac, alpha_brake]
%   alpha_trac:  acceleration scaling factor [0.5-1.0] (default 1.0)
%   alpha_brake: braking scaling factor [0.5-1.0] (default 1.0)
%
% Returns:
%   running_inter: running time for each inter-station (seconds)
%   Total_E: total energy consumption (kWh)
%   s_out, v_out, vlim_out: diagnostic output arrays

global var;           % decision variable indices per section
global vel_profile;   % speed profile (distance vs speed limit) from route
global station_info;  % station location and dwell time data
global gradient;      % gradient (%) along route
global terminal;      % terminal station info
%% Rolling stock parameters
% global Mass; % Mass of the rolling stock (kg)
global lambda;          % rotational mass coefficient
global inertial_mass;   % inertial mass equivalent (kg)
global Davis;           % Davis running resistance coefficients

%% Traction and braking characteristics (railway company-specific)
global V1_traction;     % speed transition 1 for traction (m/s): constant force until V1
global V2_traction;     % speed transition 2 for traction (m/s): constant power region
% global V3_traction;   % speed transition 3 for traction (m/s): power decreases as 1/V

global V1_brake;        % speed transition 1 for braking (m/s): constant effort until V1
global V2_brake;        % speed transition 2 for braking (m/s): constant power region
% global V3_brake;      % speed transition 3 for braking (m/s): power decreases as 1/V

global Max_tractive_power;  % maximum tractive power (kW): 3144 kW
global Max_brake_power;     % maximum braking power (W) 
global max_speed;
%global co_fric;
%global gravity;
%global max_traction;%KN
%global max_brake; %KN
%MAXIMUM ACCELERATION
global max_accel_trac;%m/s2
global max_accel_brake;%m/s2

%% ===== Adaptive acceleration/braking scaling (optional END parameters) =====
% If X has 2 extra scalars: [...base_decision_vars..., alpha_trac, alpha_brake]
% alpha values range [0.5, 1.0]. Otherwise defaults to 1.0
alpha_trac = 1.0;
alpha_brake = 1.0;
try
    nXexp = sum(var);
    if numel(X) == nXexp + 2
        alpha_trac  = min(max(X(end-1), 0.8), 1.0);
        alpha_brake = min(max(X(end  ), 0.8), 1.0);
        X = X(1:end-2);
    end
catch
end

% var=[1,3,1,3];

%% Input driving control parameters
% dwellset = zeros(1, size(station_info,1)-2);      % dwell time (not used for continuous lines)
maxspeedset     = zeros(1, size(station_info,1)-1);      % maximum speed per inter-station (m/s)
% coasting_low  = zeros(1, size(station_info,1)-1);      % lower coasting speed threshold
accrate_fixset  = zeros(1, size(station_info,1)-1);      % maximum acceleration rate (m/s²)
brakerate_fixset= zeros(1, size(station_info,1)-1);      % maximum braking rate (m/s²)
kdset           = zeros(1, size(station_info,1)-1);      % kd parameter for control modes 2/3
train_controlset= zeros(1, size(station_info,1)-1);      % train control mode per section

%% Section-wise decision variable initialization (by speed limit segment)
% Initialize coasting and cruising speed thresholds for each speed-limit section
cst_high = zeros(1, length(vel_profile)-1);  % upper coasting threshold
cst_low1 = zeros(1, length(vel_profile)-1);  % lower coasting threshold 1
cst_low2 = zeros(1, length(vel_profile)-1);  % lower coasting threshold 2
cruising = zeros(1, length(vel_profile)-1);  % target cruising speed (CC mode)
coasting = zeros(1, length(vel_profile)-1);  % target coasting speed (CR mode)

%fprintf('call from simulation file : [%s]\n', num2str(var));

for i=1:length(var)
    idx = sum(min(var(1:(i-1)), 2));  % effective DV offset: var=3 counts as 2
    switch var(i)
        case 1
            cst_high(i) = X(idx+1)/3.6;
            cst_low1(i) = cst_high(i);
            cst_low2(i) = cst_high(i);
        case 2
            cruising(i) = X(idx+1)/3.6;
            coasting(i) = X(idx+2)/3.6;
        case 3
            cst_high(i)  = X(idx+1)/3.6;
            cst_low1(i)  = X(idx+2)/3.6;
            cst_low2(i)  = cst_low1(i);  % fallback for backward sim; MID overridden adaptively in forward sim
    end
end
% % up-direction
% dwellset=[20,20,15,20,30,20,25,35,20,15,20,15,20,20];   %Dwell time
% down-direction
% dwellset=[20,20,15,20,30,20,25,35,20,15,20,15,20,20];   %Dwell time

%dwellset=dwellset(1:size(station_info,1)-2);
%dwellset(:)=30;
maxspeedset(:)=max_speed;
% coasting_low(:)=maxspeedset;
% coasting_low(:)=20;
% % maxspeedset(:)=121;
% Coasting_velset=maxspeedset;
% accrate_fixset(:)=1.2;
% accrate_fixset(:)=0.95;
accrate_fixset(:)=alpha_trac*max_accel_trac;
% brakerate_fixset(:)=2.0;
% brakerate_fixset(:)=0.57;
brakerate_fixset(:)=alpha_brake*max_accel_brake;
kdset(:)=0.5;
train_controlset(:)=1;
% Coasting_velset(:)=39/3.6;
brake_reg_efficiency=0.6;%The braking regenration efficiency is 60%
Coasting=ones(1,size(station_info,1)-1);%set a coasting flag for all inter-stations,'1'means coasting,'0' means no coasting



%% run MAIN_SIMULATION
% MAIN_SIMULATION;

% optimise the traction power calculation 2015.03.25
% large down hill is not suitable for this coasting speed motion simulator
% a coasting point motion simulator is more accurate

%2015.03.26 it is stuatable for a steep down hill, but not a perfect solution


%COMPUTE TRACTION CURVES:
% tract

%students to modify code below
%Coasting = 1;
%Coasting_vel=40/3.6;       %unit m/s
%

%DEFINE TRACK DOMAIN S is the distance variable
S_min=0;
S_max=max(vel_profile(:,1))*1000;


del_S=1;
s=S_min:del_S:S_max;
Size=length(s);
real_limit_plot = zeros(1, Size);

%DEFINE ARRAYS
vel_limit=zeros(1, Size);
vel_error=zeros(1, Size);
grad_prof=zeros(1, Size);
TractionF=zeros(1, Size);
TractionB=zeros(1, Size);
EnergyF=zeros(1, Size);
EnergyB=zeros(1, Size);
Notch=zeros(1, Size);
acceleration=zeros(1, Size);
acceler=zeros(1, Size);
accF=zeros(1, Size);
accB=zeros(1, Size);
T=zeros(1, Size);
vel=zeros(1, Size);
del_T=zeros(1, Size);
running_inter=zeros(1,length(station_info(:,1))-1);%running time of every inter-station

%DETERMINE VEL LIMIT FROM INPUT FILE
pos=2;
for i=1:Size-1
        %% %%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for STSi=1:1:(size(station_info,1)-1)
        if i*del_S/1000>=station_info(STSi,1)&& i*del_S/1000<station_info(STSi+1,1)  % i is distance, in m, but with 10m interval
            break;
        end
    end
    maxspeed=maxspeedset(STSi)/3.6;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set the speed limit according to the cruising speed of every speed
    %limit sections
    if cruising(pos-1)~=0%only for cruising sections
        vel_limit(1,i)=min(vel_profile(pos-1,2)/3.6,cruising(pos-1));  %changed to m/s
    else
        vel_limit(1,i)=vel_profile(pos-1,2)/3.6;
    end

    real_limit_plot(1,i) = vel_profile(pos-1,2)/3.6;
    
     %% %%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if maxspeed<vel_limit(1,i)
        vel_limit(1,i)=maxspeed;
    end   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    if vel_limit(1,i)>=max_speed/3.6
        vel_limit(1,i)=max_speed/3.6;
    end
    
    if i*del_S>=fix(vel_profile(pos,1)*1000)
        pos=pos+1;
    end
end

%2021.12.31 new
%DETERMINE STATION STOPS FROM INPUT FILE
if size(station_info,1)>2  %if there are other stations other than the fist and last station
    for STSi=2:(size(station_info,1)-1)
        i=fix(1000/del_S*station_info(STSi,1));
        dwell=dwellset(STSi-1);  %note this is for dwell time only
        vel_limit(1,i)=1/dwell*del_S;
        vel_limit(1,i+1)=1/dwell*del_S;
    end 
end

%DETERMINE Terminal STATION STOPS FROM INPUT FILE
pos=1;
for i=1:fix(1000/del_S*max(terminal(:,1)))
    if i*del_S>=fix(terminal(pos,1)*1000)
        vel_limit(1,i)=3/(Terminal_time)*del_S;
        vel_limit(1,i+1)=3/(Terminal_time)*del_S;
        pos=pos+1;
        %fprintf('Terminal Time : %.0f \n', Terminal_time);
    end
end

%display(terminal(pos,1));

%DETERMINE GRADIENT PROFILE FROM INPUT FILE
pos=2;
for i=1:fix(1000/del_S*max(gradient(:,1)))
    grad_prof(1,i)=-gradient(pos-1,2)/1000*9.8/(1+lambda);
    if i*del_S>=fix(gradient(pos,1)*1000)
        pos=pos+1;
    end
end

%GRADIENT SMOOTHING - not necessary
for i=5:fix(1000/del_S*max(gradient(:,1)))-4
    grad_prof(1,i)=mean(grad_prof(1,i-3:i+3));
end

%FORWARD VELOCITY CALCULATION
velF=zeros(1,Size)+0.01;

% Precompute gradient-adaptive MID per section (avoids interp1 inside hot loop)
mid_adaptive = cst_low2;
for si = 1:length(var)
    if var(si) == 3
        sec_mid_km = (vel_profile(si,1) + vel_profile(si+1,1)) / 2;
        local_g = interp1(gradient(:,1), gradient(:,2), sec_mid_km, 'linear', 'extrap');
        steepness = min(max(0, -local_g) / 10, 1.0);
        H = cst_high(si); L = cst_low1(si);
        mid_adaptive(si) = max(H - steepness * (H - L) * 0.5, L);
    end
end

Traction_power_forward=zeros(1,Size);
coast_start_flag=0;
%%
for i=1:(Size-1)
    Resistance=(Davis(1)+Davis(2)*velF(1,i)+Davis(3)*(velF(1,i))^2)/inertial_mass;
    %% %%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for STSi=1:1:(size(station_info,1)-1)
        if i*del_S/1000>=station_info(STSi,1)&& i*del_S/1000<station_info(STSi+1,1)  % i is distance, in m, but with 10m interval            break;
            break;
        end
    end
    train_control=train_controlset(STSi);
    kd=kdset(STSi);
    accrate_fix=accrate_fixset(STSi); %final speed accerlate
    vel_error(1,i)=1;
    %%%%%%%%%%%%%%%%%%%%2022.10.24 new two combine coasting-remtoring and cruising%%%%%%%%%%%%%%%
    for SPDLi=1:1:(length(vel_profile)-1)%find the current speed limit section the train
        if i*del_S/1000>=vel_profile(SPDLi,1)&& i*del_S/1000<vel_profile(SPDLi+1,1)  % i is distance, in m, but with 10m interval            break;
            %THE following value setting is to avoid change the initial
            %coast bounds changing
            LOW1=cst_low1(SPDLi);%the first lower bound of coasting-remotoring at this speed limit section
            HIGH=cst_high(SPDLi);%the higher bound of coasting-remotoring at this speed limit section
            MID=mid_adaptive(SPDLi);% precomputed: gradient-adaptive for var=3, cst_low2 otherwise
            break;
        end
    end
    %for coasting-remotoring
    if var(SPDLi)==3||var(SPDLi)==1 %coasting-remotoring
        %20230214, update, to avoid both High and low1 are higher than the
        %speed limit
        if HIGH<=LOW1%avoid the lower bound is higher than the higher bound,or both higher than speed limit
            if LOW1>=vel_limit(1,i)
                LOW1=vel_limit(1,i);
                HIGH=vel_limit(1,i);
            else
                LOW1=HIGH;
            end
        end
        % ----- NEW: clamp MID to be between [LOW1, HIGH] -----
        if exist('MID','var')
            % ensure ordering LOW1 <= MID <= HIGH
            MID = min(max(MID, LOW1), HIGH);
            % also avoid MID above speed limit
            MID = min(MID, vel_limit(1,i));
        else
            MID = HIGH;
        end

if (velF(1,i)<vel_limit(1,i+1) && velF(1,i)<MID && coast_start_flag==0)||(velF(1,i)<LOW1) %traction
                
            Traction = traction_rate(velF(1,i),max_accel_trac(1),Max_tractive_power,inertial_mass,vel_error(1,i),V1_traction,V2_traction);
            accF(1,i+1)=Traction - Resistance + grad_prof(1,i);
            
            if accF(1,i+1)>accrate_fix %determine if the acceleration is larger than the maximum value
                accF(1,i+1)=accrate_fix+(accF(1,i+1)-accrate_fix)*0.1;
                Traction=accF(1,i+1)+Resistance - grad_prof(1,i);
            end
            velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;
            %set coast_start_flag=0,this is for the situation coast_low of the next speed limit section
            %increase to a value higher than the current speed, and the current situation is
            %coasting(because the current coasting_low is lower)
            coast_start_flag=0;     
            if velF(1,i+1)>=vel_limit(1,i+1)&& HIGH~=LOW1%start coasting from the next vel_limit
                velF(1,i+1)=vel_limit(1,i+1);
                %              velF(1,i+1)=coasting_high(1,i+1);
                coast_start_flag=1;
            end
            if velF(1,i+1)>=LOW1&& HIGH==LOW1%the next step of acceleration will cause vel>coasting speed
                velF(1,i+1)=LOW1;
                
            end
            if velF(1,i+1)>=MID && HIGH~=LOW1  %start EARLY coasting from MID
                coast_start_flag=1;
            end

            % also allow coasting if reaching HIGH (backward compatibility)
            if velF(1,i+1)>=HIGH && HIGH~=LOW1
                coast_start_flag=1;
            end
        else
            if velF(1,i)==LOW1 && HIGH==LOW1%2022.02.25 new,this will cause cruising at speed coasting_low(1,i+1)/coasting_high(1,i+1)
                Traction = Resistance - grad_prof(1,i);
                %add 2022.05.30,to deal with the needed cruising traction force is
                %higer than the maximum traction.
                if Traction<=traction_rate(velF(1,i),max_accel_trac(1),Max_tractive_power,inertial_mass,vel_error(1,i),V1_traction,V2_traction)
                    accF(1,i+1)=0;
                    %only when the speed limit is smaller than coasting speed,set vel=speed limit
                    velF(1,i+1)=min(vel_limit(1,i+1),LOW1);
                else
                    Traction=traction_rate(velF(1,i),max_accel_trac(1),Max_tractive_power,inertial_mass,vel_error(1,i),V1_traction,V2_traction);
                    accF(1,i+1)=Traction - Resistance + grad_prof(1,i);
                    velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                coast_start_flag=0;
            else    %coasting_high ~=coasting_low
                %normal coasting, after traction to the speed limit or
                %coasting_high
                if velF(1,i)<vel_limit(1,i+1) && coast_start_flag==1
                    
                    Traction = 0;
                    accF(1,i+1)= - Resistance + grad_prof(1,i);
                    velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;%add 20220504,must be assigned before use
                    %encounter steep slope when coasting,use brake force to
                    %keep to the speed limit
                    if velF(1,i+1)>=vel_limit(1,i+1)
                        Traction =  Resistance - grad_prof(1,i);
                        accF(1,i+1)=0;
                    end
                    velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;
                    %%%coasting to the lower coasting speed, then stop coasting
                    if velF(1,i+1)<=LOW1
                        coast_start_flag=0;
                    end
                else
                    if vel_limit(1,i+1)<1 %arrive the station,cruising untile depart
                        Traction= Resistance - grad_prof(1,i);
                        accF(1,i+1)=0;
                        velF(1,i+1)=vel_limit(1,i+1);
                        coast_start_flag=0; %set coasting flag=0 after the station
                        %                   end
                    else
                        % This including two situations caused by speed limit decrease
                        %one is:after decrease,v(i)>v_limit(i+1)
                        %Another situation : when train coasting to coasting_low,speed
                        %limit decrease at this position. the speed is lower
                        %than speed limit of the next section and higher than
                        % coasting_high of the next section.then we let the
                        %train continue coasting to coasting_low of the next
                        %speed limit section
                        Traction = 0;
                        accF(1,i+1)= - Resistance + grad_prof(1,i);
                        %when v(i)>v_limit(i+1), encounter steep slope when coasting,use brake force to
                        %keep to the speed limit
                        if accF(1,i+1)>0
                            Traction =  Resistance - grad_prof(1,i);%20220512 revised, the former one the wrong as ' Traction = -Resistance + grad_prof(1,i)'
                            accF(1,i+1)=0;
                        end
                        velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;
                        if velF(1,i+1)>vel_limit(1,i+1)
                            velF(1,i+1)=vel_limit(1,i+1);
                        end
                        if  coast_start_flag==0
                            coast_start_flag=1; %when the steep slope disappear, continue coasting
                        end
                        
                    end
                end
            end
        end
    else%for cruising var=2
        %%%%%%%%%%%%2021.12.15 new,calculate traction according to different stages%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Traction = traction_rate(velF(1,i),max_accel_trac(1),Max_tractive_power,inertial_mass,vel_error(1,i),V1_traction,V2_traction);
        %%%%%%%%%%%%2021.12.15 new,from the following 'if-else'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        accF(1,i+1)=Traction - Resistance + grad_prof(1,i);
        % NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if accF(1,i+1)>accrate_fix                                               %NEW
            accF(1,i+1)=accrate_fix+(accF(1,i+1)-accrate_fix)*0.1;                   %NEW
            Traction=accF(1,i+1)+Resistance - grad_prof(1,i);
        end                                                                   %NEW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if velF(1,i)<vel_limit(1,i+1)
            velF(1,i+1)=(velF(1,i).^2+2*accF(1,i+1)*(del_S)).^0.5;
            if velF(1,i+1)>vel_limit(1,i+1)
                velF(1,i+1)=vel_limit(1,i+1);
            end
            Traction_power_forward(1,i)=Traction*inertial_mass*1000*velF(1,i);
            
        else
            %This if condition is to allow for deceleration at full power
            %while assending a gradient
            if (velF(1,i)==vel_limit(1,i+1))&&(accF(1,i+1)<0)  %added in 2021.12.17,to deal with speed restriction decrease
                %         if accF(1,i+1)<0    %delete in 2021.12.17
                velF(1,i+1) = (velF(1,i)^2+2*(accF(1,i+1))*(del_S))^0.5;
                Traction_power_forward(1,i)=Traction*inertial_mass*1000*velF(1,i);
            else    %%%%%%cruising
                accF(1,i+1)=0;
                velF(1,i+1) = vel_limit(1,i+1);
                Traction=Resistance-grad_prof(1,i);
                Traction_power_forward(1,i)=Traction*inertial_mass*1000*velF(1,i);
            end
            
        end
        
    end
    Traction_power_forward(1,i)=Traction*inertial_mass*1000*velF(1,i);
%     if ~isreal(velF(1,i))
%         a=1;
%     end
end

%%
Traction_power_backward=zeros(1,Size);
%%

%BACKWARD VELOCITY CALCULATION
cst_flag=zeros(1,length(vel_profile)-1);%the coasting flag of every speed limit section
velB=zeros(1,Size)+0.01;
for i=(Size-1):-1:1
    if velB(1,i+1)<vel_limit(1,i)
        
        %% %%%%% NEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for STSi=1:1:(size(station_info,1)-1)
            if i*del_S/1000>=station_info(STSi,1)&& i*del_S/1000<station_info(STSi+1,1)
                break;
            end
        end
        brakerate_fix=brakerate_fixset(STSi);
        for SPDLi=1:1:(length(vel_profile)-1)%find the current speed limit section the train
            if i*del_S>vel_profile(SPDLi,1)*1000&& i*del_S<=vel_profile(SPDLi+1,1)*1000  % i is distance, in m, but with 10m interval            break;
                break;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Resistance=(Davis(1)+Davis(2)*velB(1,i+1)+Davis(3)*(velB(1,i+1))^2)/inertial_mass;
        if var(SPDLi)==1%backward braking to the speed limit
            brake_rate=brake(velB(1,i+1), max_accel_brake, Max_brake_power, inertial_mass,V1_brake,V2_brake);%, Coasting, Coasting_vel);
            accB(1,i)=Resistance + brake_rate- grad_prof(1,i+1);
        else  %cruising coasting mode and coasting-remotoring with two lower bounds
            if var(SPDLi)==2
                Coasting_vel=coasting(SPDLi);
            else
                Coasting_vel=cst_low2(SPDLi);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            brake_rate=brake(velB(1,i+1), max_accel_brake, Max_brake_power, inertial_mass,V1_brake,V2_brake);%, Coasting, Coasting_vel);
            %% new2022.10.22 multi coasting speed
            %the cst_flag is 0 at first, when the speed exceed the coasting speed for the first time,
            %this flag is set to 1.Then the train will keep coasting regardless
            %of whether the speed is smaller than the coasting speed. This is
            %to designed to deal with the speed fluctuation under the coasting speed in steep downhill
            %but this will cause the speed coasting to 0 when the steep
            %downhill is too long. this problem haven't been solved here
            if cst_flag(SPDLi)==0 && velB(1,i+1)>Coasting_vel
                cst_flag(SPDLi)=1;
                brake_rate=0;
            end
            if cst_flag(SPDLi)==1
                %to avoid speed reduce to 0, let the train cruising when the
                %speed is less than 20km/h
                if velB(1,i+1)>=20/3.6
                    brake_rate=0;
                else
                    brake_rate=grad_prof(1,i+1)-Resistance;
                end
            end
            accB(1,i)=Resistance + brake_rate- grad_prof(1,i+1);
        end
        if accB(1,i)>brakerate_fix
            accB(1,i)=brakerate_fix+(accB(1,i)-brakerate_fix)*0.1;
            brake_rate=accB(1,i)-Resistance +grad_prof(1,i+1);
        end
        Traction_power_backward(1,i+1)=brake_rate*inertial_mass*1000*velB(1,i+1);
        velB(1,i)=(velB(1,i+1)^2+2*accB(1,i)*(del_S))^0.5;
        % velB reach vel_limit, velB cannot become vel,because of velF can
        if velB(1,i)>=vel_limit(1,i)
            velB(1,i)=vel_limit(1,i);
        end
        
    else
        velB(1,i) = vel_limit(1,i);
    end
    if ~isreal(velB(1,i))
        a=1;
    end
end


TractionF_new=zeros(1,Size);
TractionB_new=zeros(1,Size);
%ACTUAL VELOCITY PROFILE %%  select the lower speed
for i=1:Size
    if velF(1,i)<=velB(1,i)
        acceler(1,i)=accF(1,i);
        vel(1,i)=velF(1,i);
        TractionF_new(1,i)=Traction_power_forward(1,i);
        
        if TractionF_new(1,i)<0                          %NEW sometimes it's braking when cruising
            TractionB_new(1,i)=-TractionF_new(1,i);       %NEW
            TractionF_new(1,i)=0;                        %NEW
        end
        
    else
        acceler(1,i)=-accB(1,i);
        vel(1,i)=velB(1,i);
        TractionB_new(1,i)=Traction_power_backward(1,i);
        
        %% 2022.01.03 new,when there is steep slope during coasting, pulling is needed to change to cruising
        if TractionB_new(1,i)<0
            TractionF_new(1,i)=-TractionB_new(1,i);
            TractionB_new(1,i)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%TIME CALCULATION
for i=1:Size-1
    del_T(1,i+1)=2*del_S/(vel(1,i)+vel(1,i+1));
    if vel(1,i)==0
        T(1,i+1)=T(1,i)+1;
    else
        T(1,i+1)=T(1,i)+del_T(1,i+1);
    end
end

%% 2022.02.18 calculate the number of stations
% to get rid of the '0'station position value which aims to
% fill in the blank in the end of the station info
if station_info(1,1)==0
    station_num=length(nonzeros(station_info(:,1)))+1;
else
    station_num=length(nonzeros(station_info(:,1)));
end

%2021.12.28 inter-station TIME CALCULATION
%undated 2022.02.18 
%get ride of the journey time does not match with the totoal running time plusing the total dwell time 
for STSi=2:1:station_num %(STSi-1)is the inter-station
    if STSi>2 && STSi<station_num
        %this is for inter-stations,not the fist inter-station and not the
        %last inter-station
       %inter-station running time=T(Train reaches the STSi station)-T(Train leaves the STSi-1 station)
       running_inter(STSi-1)=T(1,fix(station_info(STSi,1)*1000/del_S))-T(1,fix(station_info(STSi-1,1)*1000/del_S)+1);
   else
       if STSi==2 && station_num~=2 
           %this is for the first inter-stationand the jouney has more than one inter-station
           running_inter(STSi-1)=T(1,fix(station_info(STSi,1)*1000/del_S)); 
       else
           if STSi==station_num && station_num~=2 
           %for the last inter-station and the jouney has more than
           %one inter-station
            running_inter(STSi-1)=T(1,fix(station_info(STSi,1)*1000/del_S)+1)-T(1,fix(station_info(STSi-1,1)*1000/del_S)+1); 
           else
               %this is for the first inter-stationand the jouney has only one inter-station
               running_inter(STSi-1)=T(1,fix(station_info(STSi,1)*1000/del_S)+1); 
           end
       end
    end
end

%TRACTION (F-forward B-backward) CALCULATION
for i=2:Size-2
    Resistance=(Davis(1)+Davis(2)*vel(1,i)+Davis(3)*(vel(1,i))^2)/inertial_mass;
    Notch(1,i)=((acceler(1,i-1)+acceler(1,i-1))/2+ Resistance - grad_prof(1,i))*vel(1,i)*inertial_mass*1000;
    if Notch(1,i)>0
        TractionB(1,i)=0;
    else
        TractionB(1,i)=-Notch(1,i);
    end
    Notch(1,i)=((acceler(1,i+1)+acceler(1,i+1))/2 + Resistance - grad_prof(1,i))*vel(1,i)*inertial_mass*1000;
    if Notch(1,i)>0
        TractionF(1,i)=Notch(1,i);
    else
        TractionF(1,i)=0;
    end
    del_T(1,i)=2*del_S/(vel(1,i)+vel(1,i+1));
end

EnergyF_new=zeros(1,Size);
EnergyB_new=zeros(1,Size);



%Energy (F-forward B-backward) CALCULATION
for i=2:Size-2
    EnergyF(1,i)=TractionF(1,i)/vel(1,i)*del_S;
    EnergyF_new(1,i)=TractionF_new(1,i)/vel(1,i)*del_S;
    if TractionB(1,i)> Max_brake_power
        EnergyB(1,i)=Max_brake_power/vel(1,i)*del_S;    
    else
        EnergyB(1,i)=TractionB(1,i)/vel(1,i)*del_S;
    end
    EnergyB_new(1,i)=TractionB_new(1,i)/vel(1,i)*del_S;
end

% runs plot files
%plotter_all
altitude_calculate;  %altitude calculation
%kWh_Mech_F=sum(EnergyF)/3.6/1000000;
%kWh_Mech_B=sum(EnergyB)/3.6/1000000;
%journey_time=max(T);


%kWh_Mech_F_new=sum(EnergyF_new)/3.6/1000000;
%kWh_Mech_B_new=sum(EnergyB_new)/3.6/1000000;


%energy_consumed=(cumsum(EnergyF)-brake_reg_efficiency*cumsum(EnergyB))/3.6/1000000;%2021.12.21 add regeneration efficiency
%kinetic_energy=1/2*inertial_mass.*vel.*vel*1000/3.6/1000000;
%potential_energy=altitude*9.81*(Mass)/3.6/1000;
%resistance_energy= cumsum((Davis(1)+Davis(2).*vel+Davis(3).*(vel.^2)).*del_S)*1000/3.6/1000000;

Total_E=sum(EnergyF_new)/1000/3600;

% --- [PATCH] OUTPUT UNTUK GRAFIK ---
    s_out = s;              % Data Jarak
    v_out = vel;            % Data Kecepatan Aktual
    %vlim_out = vel_limit;   % Data Speed Limit
    vlim_out = real_limit_plot;
    % -----------------------------------

% figure()
% plot(s,velF);
% hold on
% plot(s,velB);
% figure()
% plot(T,vel);
% figure()
% plot(s/1000,vel_limit*3.6)
% figure()
% plot(s/1000,velB*3.6)
% hold on
% plot(s/1000,velF*3.6)
% hold on
% plot(s/1000,vel_limit*3.6,'--')
% hold on
% % plot(s/1000,cst_high*ones(1,length(s)),'--')
% % hold on
% % plot(s/1000,cst_low*ones(1,length(s)),'--')
% xlabel('distance (km)')
% ylabel('speed (km/h)')
% legend('SpeedB', 'SpeedF','speed limit','upper bound of coasting','lower bound of coasting')
% A=[Total_E running_inter X];
% 
% if Total_E~=100000
%     writematrix(A,'energy_time_GA_CR.xls','WriteMode','append');
% end
%% This is for the surrogateopt that don't have constrain function
% if running_inter<199 || running_inter>201
%     return
% end
% A=[Total_E running_inter X];
% if Total_E~=1000
%     writematrix(A,'energy_time_GA.xls','WriteMode','append');
% end

end

function [FA] = traction_rate(vel, max_accel, Power, Mass,vel_error,V1_traction,V2_traction)%here Mass is the effective mass
    %2021.12.15 new,calculate traction according to different stages
    if vel<=V1_traction
        FA=Power/V1_traction*vel_error/(Mass*1000);
    else
        if (vel>V1_traction)&&(vel<=V2_traction)
            FA=Power*vel_error/vel/(Mass*1000);
        else
            FA=Power*V2_traction*vel_error/vel^2/(Mass*1000);
        end
    end
    if FA>max_accel %%%make sure the train doesn't exceed the maximum acceleration
        FA=max_accel;
    end
end

function [b] = brake(vel, max_accel, Power, Mass,V1_brake, V2_brake)%, Coasting, Coasting_vel),Input parameter ('Mass') of this function can be the effective mass.
    if vel>1 
        %%%%%%%%%%%%%%%%%%%2021.12.16 new, brake with different stages%%%%%%%%%
        if vel<=V1_brake
            b=Power/V1_brake/(Mass*1000);
        else
            if (vel>V1_brake)&&(vel<=V2_brake)
                b=Power/vel/(Mass*1000);
            else
                b=Power*V2_brake/vel^2/(Mass*1000);
            end
        end
    %     b=Power/Mass/1000/vel;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if b>max_accel;
            b=max_accel;
        end
    else
        b=max_accel;
    end
end

