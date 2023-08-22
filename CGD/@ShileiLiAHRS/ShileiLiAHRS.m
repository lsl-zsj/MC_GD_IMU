classdef ShileiLiAHRS < handle
%ShileiLiAHRS Implementation of ShileiLi's IMU and AHRS algorithms
%
%   For more information see: 
%
%   Date          Author          Notes
%   28/09/2011    ShileiLi    Initial release

%% Public properties
    properties (Access = public)
        cnt=0;
        SamplePeriod = 1/256;
        Quaternion = [1 0 0 0];     % output quaternion describing the Earth relative to the sensor
        Beta = 1;               	% algorithm gain
        Sigma_g=1;
        Sigma_m=1;
        Err=zeros(6,1);
    end

    %% Public methods
    methods (Access = public)
        function obj = ShileiLiAHRS(varargin)
            obj.cnt=0;
            for i = 1:2:nargin
                if  strcmp(varargin{i}, 'SamplePeriod'), obj.SamplePeriod = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Quaternion'), obj.Quaternion = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Beta'), obj.Beta = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Sigma_g'), obj.Sigma_g = varargin{i+1};
                elseif  strcmp(varargin{i}, 'Sigma_m'), obj.Sigma_m = varargin{i+1};
                else error('Invalid argument');
                end
            end
        end
        function obj = Update(obj, Gyroscope, Accelerometer, Magnetometer)
            obj.cnt=obj.cnt+1;
            % obtain the initial oriention using ecompass algorithm
            if(obj.cnt==1)
               accr=Accelerometer;
               gyror=Gyroscope;
               magr=Magnetometer;
               % NED coordinate
               r_down=accr';
               r_east=cross(accr',magr');
               r_north= cross(r_east, r_down);

               r_down=r_down/norm(r_down);
               r_east=r_east/norm(r_east);
               r_north=r_north/norm(r_north);

               R_=[r_north,r_east,r_down];
               Q__ = quaternion(R_, 'rotmat', 'frame');    
               obj.Quaternion=compact(Q__);
            end
            q = obj.Quaternion; % short name local variable for readability

            % Normalise accelerometer measurement
            if(norm(Accelerometer) == 0), return; end	% handle NaN
            Accelerometer = Accelerometer / norm(Accelerometer);	% normalise magnitude

            % Normalise magnetometer measurement
            if(norm(Magnetometer) == 0), return; end	% handle NaN
            Magnetometer = Magnetometer / norm(Magnetometer);	% normalise magnitude

            % Reference direction of Earth's magnetic feild
            h = quaternProd(q, quaternProd([0 Magnetometer], quaternConj(q)));
            b = [0 norm([h(2) h(3)]) 0 h(4)];

            % Gradient decent algorithm corrective step
            F = [2*(q(2)*q(4) - q(1)*q(3)) - Accelerometer(1)
                2*(q(1)*q(2) + q(3)*q(4)) - Accelerometer(2)
                2*(0.5 - q(2)^2 - q(3)^2) - Accelerometer(3)
                2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3)) - Magnetometer(1)
                2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4)) - Magnetometer(2)
                2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2) - Magnetometer(3)];
            J = [-2*q(3),                 	2*q(4),                    -2*q(1),                         2*q(2)
                2*q(2),                 	2*q(1),                    	2*q(4),                         2*q(3)
                0,                         -4*q(2),                    -4*q(3),                         0
                -2*b(4)*q(3),               2*b(4)*q(4),               -4*b(2)*q(3)-2*b(4)*q(1),       -4*b(2)*q(4)+2*b(4)*q(2)
                -2*b(2)*q(4)+2*b(4)*q(2),	2*b(2)*q(3)+2*b(4)*q(1),	2*b(2)*q(2)+2*b(4)*q(4),       -2*b(2)*q(1)+2*b(4)*q(3)
                2*b(2)*q(3),                2*b(2)*q(4)-4*b(4)*q(2),	2*b(2)*q(1)-4*b(4)*q(3),        2*b(2)*q(2)];
            
            error=F;
            % correntropy
            GF=ones(size(F));
            % accelerometer
            for i=1:3
               GF(i)= exp(-(F(i)^2)/(2*obj.Sigma_g^2));
            end
            % magnerometer
            for i=4:6
               GF(i)= exp(-(F(i)^2)/(2*obj.Sigma_m^2));
            end
            
            % step old
            stepori = J'*F;
            step_norm = norm(stepori);	% normalise step magnitude
            % step new
            F=F.*GF; % shrinkage
            step=J'*F;
            step=step/step_norm;
            
            % Compute rate of change of quaternion
            qDot = 0.5 * quaternProd(q, [0 Gyroscope(1) Gyroscope(2) Gyroscope(3)]) - obj.Beta * step';

            % Integrate to yield quaternion
            q = q + qDot * obj.SamplePeriod;
            
            if(q(1)<0)
            q=-q;
            end
            obj.Quaternion = q / norm(q); % normalise quaternion
            obj.Err=error;
        end
        function obj = UpdateIMU(obj, Gyroscope, Accelerometer)
            q = obj.Quaternion; % short name local variable for readability

            % Normalise accelerometer measurement
            if(norm(Accelerometer) == 0), return; end	% handle NaN
            Accelerometer = Accelerometer / norm(Accelerometer);	% normalise magnitude

            % Gradient decent algorithm corrective step
            F = [2*(q(2)*q(4) - q(1)*q(3)) - Accelerometer(1)
                2*(q(1)*q(2) + q(3)*q(4)) - Accelerometer(2)
                2*(0.5 - q(2)^2 - q(3)^2) - Accelerometer(3)];
            J = [-2*q(3),	2*q(4),    -2*q(1),	2*q(2)
                2*q(2),     2*q(1),     2*q(4),	2*q(3)
                0,         -4*q(2),    -4*q(3),	0    ];
            step = (J'*F);
            step = step / norm(step);	% normalise step magnitude

            % Compute rate of change of quaternion
            qDot = 0.5 * quaternProd(q, [0 Gyroscope(1) Gyroscope(2) Gyroscope(3)]) - obj.Beta * step';

            % Integrate to yield quaternion
            q = q + qDot * obj.SamplePeriod;
            obj.Quaternion = q / norm(q); % normalise quaternion
        end
    end
end