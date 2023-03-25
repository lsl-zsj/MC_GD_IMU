function out =EMDI_C(Acc,Gyr,Mag,sample_freq, TauAcc, TauMag, Zeta, accRating,Sigma_acc,Sigma_mag)

% #############################################################################################################################################################
% IMU ORIENTATION ESTIMATION
% T. Seel, S. Ruppin. Eliminating the Effect of Magnetic Disturbances on the Inclination Estimates of Inertial Sensors. IFAC-PapersOnLine,
% Contact: {seel,ruppin}@control.tu-berlin.de
% 
% All units are SI units.
%    Inputs:
%    estimation_state:  state of the estimation before the update
%                       size of estimation_state [1 18]
%                       estimation_state = [quaternion, estimatedBias, validAccDataCount, ratingWindow]
%                       length(quaternion)          == 4
%                       length(estimatedBias)       == 4
%                       length(validAccDataCount)   == 1
%                       length(ratingWindow)        == windowLength
%   Acc:                measurement of the accelerometer [x,y,z][m/s^2]
%                       size == [1 3]
%   Gyro:               measurement of the gyro          [x,y,z][rad/s]
%                       size == [1 3]
%   Mag:                measurement of the magnetometer  [x,y,z][any unit] 
%                       size == [1 3]
%   sample_freq:        sample frequency of measurement data
%   TauAcc, TauMag:     time constants for correction (50% time) [must be >0]
%                       [in seconds]                       
%   Zeta:               bias estimation strength [no unit]
%                       0: do not estimate gyroscope bias
%                       >0: estimate gyro bias (larger means more aggressive)
%                       Recommended range: 0–5
%   useMeasRating:      reduce weight/trust of Acc (i.e. increase TauAcc) when Acc norm is not close to 9.81m/s^2
%                       0: do not use raw rating of accelerometer
%                       >0: use raw rating of acclerometer (larger means more aggressive)
%                       Recommended range: 1–10
%   Outputs:
%    q?                quaternion
%   bias:               bias 
%   eIncl:              error angle of inclination [scalar]
%   eAzi:               error angle of azimuth [scalar]
% #############################################################################################################################################################

% #### Init declarations #####
  len=length(Acc);
  bias=0;
  windowLength    = 10;   
  window = zeros(1,windowLength);
  zeroEps                         = 1D-6;                                                            % treshold for if zero statementes

  q=zeros(len,4);
  b=zeros(len,3);
  inc=zeros(len,1);
  azi=zeros(len,1);

  


  for i=1:len
   % 
%    validAccDataCount=i-1;
%    if (validAccDataCount+1)/sample_freq/2 < TauAcc/2                                                   % start bias estimation after TauAcc/2 is reached
%        Zeta = 0;                                                                                       % to prevent estimating the bias out of the error due 
%    else
%        Zeta = Zeta1;
%    end 
%    TauAcc = min((validAccDataCount)/sample_freq/2,TauAcc1);                                           % Speed-up acceleromter-based correction in the first samples
%    TauMag = min((validAccDataCount)/sample_freq/2,TauMag1); 
    
   % obtain the current data
   acc=Acc(i,:);
   gyr=Gyr(i,:);
   mag=Mag(i,:);
   
   errorAngleIncl=0;
   errorAngleAzi=0;
   %% acc rating
   correctionRating                = 1;                                  % default:no down rating of correction gain
   if (accRating > 0)                                                    % check if rating is enabled
       window=circshift(window,-1);                                      % shift window/buffer one sample up
       window(windowLength)        = abs(norm(acc)-9.81);                % add new sample to buffer
       correctionRating            = 1/(1+max(window)*accRating);     % calculate correction rating (1: normal thrust, <1 down rating)
   end
   %% initial quaternion
   if(i==1)
       % NED coordinate
       r_down=acc';
       r_east=cross(acc',mag');
       r_north= cross(r_east, r_down);

       r_down=r_down/norm(r_down);
       r_east=r_east/norm(r_east);
       r_north=r_north/norm(r_north);

       R_=[r_north,r_east,r_down];
       Q__ = quaternion(R_, 'rotmat', 'frame');    
       q_gyro_=compact(Q__);
   end
  %% obtain the quaternion by gyroscope
    gyr                         = gyr+ bias ;
    gyro_norm                   = norm(gyr);                                                   % get norm to calculate prediction angle
    prediction_ang              = gyro_norm / (sample_freq);                                        % calculate angle to correct in angle/axis
    dq_gyro                     = [cos(prediction_ang/2), ...                                       % specify quaternion for gyro correctin
                                   sin(prediction_ang/2)*(gyr/gyro_norm)];                     % "
    q_gyro                      = quaternionMultiply(q_gyro_, dq_gyro);                              % correct initial quaternion by dq_gyro in IMU frame!
    
    %% 
    % #### Accelerometer-based correction ####
    % in: q_gyro   out: q_gyro_acc
    q_gyro_acc  = q_gyro;                                                                               % define to get an output even if acc is zero
    gravref_fixedframe              = [0 0 1];                                                          % define gravitation reference (vertical) in fixed frame
    gravref_imuframe                = quaternionCoordTransform(q_gyro,[0,gravref_fixedframe]);          % transform gravitation reference into IMU frame

    if (acc(1) ~= 0 || acc(2) ~= 0 || acc(3) ~= 0)                                          % check acc to be valid
        acc = acc/norm(acc);                                                                % normalize acc
        if (abs(acc*gravref_imuframe(2:4)')< 1)                                                     % perform correction only if acc and gravref_imuframe do NOT coincide
            errorAngleIncl          = acos(acc*gravref_imuframe(2:4)');                             % calculate error between reference and measurment          
            kp_acc                  = correctionRating ...                                              % calculate correction gain kp_acc from time constant
                                      * (1 - 1.4*TauAcc*sample_freq/(1.4*TauAcc*sample_freq+1));        % "
            
            errorAngleIncl =    errorAngleIncl*exp(-(errorAngleIncl)^2/(2*Sigma_acc^2));                  
            
            correction_ang          = kp_acc * errorAngleIncl;                                          % calculate angle for correction

            correctionaxis_imuframe = cross(acc, gravref_imuframe(2:4));                            % rotation axis of correction is perpendicular to measurement and reference
            correctionaxis_imuframe = correctionaxis_imuframe/norm(correctionaxis_imuframe);            % normalize rotation axis
            dq_acc                  = [cos(correction_ang/2), ...                                       % build correction quaternion from axis and angle
                                       sin(correction_ang/2)*correctionaxis_imuframe];                  % "
            q_gyro_acc              = quaternionMultiply(q_gyro,dq_acc);                                % correct gyro corrected quaternion by acc correction
            if kp_acc < 1
                ki_acc                  = (Zeta^2/160)*1.4*sample_freq *kp_acc^2/(1-kp_acc);            % calculate ki_acc for bias estimation
                bias                    = bias + ...                                                    % estimate bias from accelerometer correction
                                          ki_acc * errorAngleIncl*correctionaxis_imuframe;              % "
            end
            if (abs(norm(gravref_imuframe)-1)>zeroEps)                                                  % check for error (optional)
                disp('gravref_imuframe is not a unit quaternion'); end                                  % "
            if (abs(norm(dq_acc)-1)>zeroEps)                                                            % "
                disp('dq_acc is not a unit quaternion'); end                                            % "
            if (abs(norm(q_gyro_acc)-1)>zeroEps)                                                        % "
                disp('q_gyro_acc is not a unit quaternion'); end                                        % "
        end
    end
    %%
    % #############################################################################################################################################################
    % #### Magnetometer-based correction ####
    % in: q_gyro_acc   out: q_gyro_acc_mag
    q_gyro_acc_mag                  = q_gyro_acc;                                                       % define to get an output even if input is zero
    if (mag(1) ~= 0 || mag(2) ~= 0 || mag(3) ~= 0)                                          % check mag to be valid
        magref_fixedframe           = [1 0 0];                                                          % define magnetic field reference in fixed frame (choose [+/-1 0 0] or [0 +/-1 0] to define which fixed-frame coordinate axis points north/south)
        magref_imuframe             = quaternionCoordTransform(q_gyro_acc, [0,magref_fixedframe]);      % transform magnetic field reference into IMU frame
        gravref_imuframe            = quaternionCoordTransform(q_gyro_acc, [0, gravref_fixedframe]);    % recalculate IMU frame coordinates of vertical axis (might be skipped)

        if (abs(mag*gravref_imuframe(2:4)')<norm(mag))                                          % if measured magnetic field is vertical, perform NO correction
            mag_projected       = (mag' - (gravref_imuframe(2:4)*mag') ...                  % projection of magnetic measurement into horizontal plane
                                                   *gravref_imuframe(2:4)')';                           % "
            mag_projected       = mag_projected/norm(mag_projected);                        % normalize projection
            if (abs(mag_projected*magref_imuframe(2:4)')<1)                                         % if projected measurement and reference agree, perform NO correction
                errorAngleAzi       = acos(mag_projected*magref_imuframe(2:4)');                    % calculate error angle between reference and measurment
                kp_mag              = (1 - 1.4*TauMag*sample_freq/(1.4*TauMag*sample_freq+1));          % calculate correction gain  kp_mag from time constant
                
                errorAngleAzi = errorAngleAzi* exp(-(errorAngleAzi)^2/(2*Sigma_mag^2));
                correction_ang      = kp_mag*errorAngleAzi;                                            % calculate angle for correction
                correctionaxis_imuframe     = cross(mag_projected,magref_imuframe(2:4));            % rotation axis of correction is perpendicular to measurement and reference
                correctionaxis_imuframe     = correctionaxis_imuframe/norm(correctionaxis_imuframe);    % normalize rotation axis
                dq_mag              = [cos(correction_ang/2), ...                                       % build correction quaternion from axis and angle
                                       sin(correction_ang/2)*correctionaxis_imuframe];                  % "
                q_gyro_acc_mag      = quaternionMultiply(q_gyro_acc,dq_mag);                            % correct gyro&acc-corrected quaternion by mag correction
                if kp_mag < 1
                    ki_mag              = (Zeta^2/160)*1.4*sample_freq *kp_mag^2/(1-kp_mag);            % calculate ki_mag for bias estimation
                    bias                = bias + ...                                                    % estimate bias from magnetometer correction
                                          ki_mag * errorAngleAzi * correctionaxis_imuframe;             % "
                end
                                  
                if (abs(norm(dq_mag)-1)>zeroEps)                                                        % check for error (optional)
                disp('dq_mag is not a unit quaternion'); end                                            % "
            end
        end
    end
    if(q_gyro_acc_mag(1)<0)
        q_gyro_acc_mag=-q_gyro_acc_mag;
    end
    q_gyro_acc_mag= q_gyro_acc_mag/norm(q_gyro_acc_mag);                                            % if not unit quaternion, normalize it to prevent error porpagation
    q_gyro_=q_gyro_acc_mag;
    % store the data
    
    q(i,:)=q_gyro_acc_mag;
    b(i,:)=bias;
    inc(i)=errorAngleIncl;
    azi(i)=errorAngleAzi;
    
  end
 
    out.q=q;
    out.b=b;
    out.inc=inc;
    out.azi=azi;
  
end

function [q3] = quaternionMultiply (q1,q2)
    % q1 and q2 have to be Nx4 or 1x4
    q3 = [];
    q3(:,1) = q1(:,1).*q2(:,1) - q1(:,2).*q2(:,2) - q1(:,3).*q2(:,3) - q1(:,4).*q2(:,4);
    q3(:,2) = q1(:,1).*q2(:,2) + q1(:,2).*q2(:,1) + q1(:,3).*q2(:,4) - q1(:,4).*q2(:,3);
    q3(:,3) = q1(:,1).*q2(:,3) - q1(:,2).*q2(:,4) + q1(:,3).*q2(:,1) + q1(:,4).*q2(:,2);
    q3(:,4) = q1(:,1).*q2(:,4) + q1(:,2).*q2(:,3) - q1(:,3).*q2(:,2) + q1(:,4).*q2(:,1);
    
end

function [qInvert] = quaternionInvert (q)
    if (size(q, 2) ~= 4)
        error('input has to be Nx4')
    end
    qInvert = [q(:,1) -q(:,2:4)];
end
function [q3] = quaternionCoordTransform (q1,q2)
   % This function will do a rotation with two quaternions
   % Result: q1' * q2 * q1
   % The result will always be a quaternion [x x x x]
   q1Inv = quaternionInvert(q1);
   q3 = quaternionMultiply(q1Inv,q2);
   q3 = quaternionMultiply(q3,q1);
end