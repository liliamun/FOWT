%  Analysis of a Cable using Truss Element
%
clear all
clc

NE = 10;   % Number of Elements
NN = NE + 1;  
NC = NE + 1; % Number of Nodes
NDOF = 2 * NN;
EA = 1.69E5;
w = 3.2020;
LX = 704.197; LZ = 320.5;  L0 = sqrt(LX*LX+LZ*LZ);
LX1 = LX+5; L1 = sqrt(LX1*LX1+LZ*LZ);
Tref = EA*(L1-L0)/L0;
FX = Tref*LX1/L1;

% Initial Node Coordinates
for i = 1:NN
    XG(i) = (i-1)*LX1/NE;
    ZG(i) = (i-1)*LZ/NE;
end
XG0 = XG; ZG0 = ZG;

% Initial Element Length Calculation
for i=1:NE
    X1 = XG(i); X2 = XG(i+1);
    Z1 = ZG(i); Z2 = ZG(i+1);
    L_0(i) = sqrt((X2 - X1)^2 + (Z2 - Z1)^2);
end

plot(XG, ZG, 'b', 'LineWidth', 2); 
hold on;

% Initialize force vectors
FORCE_Mem = zeros(NE,1); 
iter = 1;
CN = [];

% for i = 1:NE
%     FORCE_Mem(i) = Tref;
% end

for density = 1.0 : 1.0 : 8
    error = 100.0;
        while(error > 0.01)
        
        XK = zeros(NDOF, NDOF); 
        XKG = zeros(NDOF, NDOF); 
        XFORCE = zeros(NDOF, 1); 
        DISP = zeros(NDOF, 1); 
        XFORCE2 = zeros(NDOF, 1);
        
        for i = 1:NE

            X1 = XG(i); X2 = XG(i+1);
            Z1 = ZG(i); Z2 = ZG(i+1);
            IDOF = [2*i-1; 2*i; 2*i+1; 2*i+2]; % Degrees of freedom
            LENGTH(i) = sqrt((X2 - X1)^2 + (Z2 - Z1)^2);
            dL = (X2 - X1) / LENGTH(i); 
            dN = (Z2 - Z1) / LENGTH(i);
            T = [dL dN 0 0; -dN dL 0 0; 0 0 dL dN; 0 0 -dN dL];  % Transformation
            
            if iter == 1
                FORCE_Mem(i) = Tref;
            else
                FORCE_Mem(i) = EA * (LENGTH(i) - L_0(i)) / L_0(i) + Tref;  
            end

            EF2 = FORCE_Mem(i) * T' * [-1; 0; 1; 0];
            FORCE_Mem_X(i) = FORCE_Mem(i) * dL;
            FORCE_Mem_Z(i) = FORCE_Mem(i) * dN;
            EF = -density * L_0(i) * [0 0.5 0 0.5];
            XKE1 = EA * [1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0] / LENGTH(i);
            XGE1 = FORCE_Mem(i) * [0 0 0 0; 0 1 0 -1; 0 0 0 0; 0 -1 0 1] / LENGTH(i);
            XKE = T' * XKE1 * T;      
            XKGE = T' * XGE1 * T; 

            for j = 1:4
               XFORCE(IDOF(j)) = XFORCE(IDOF(j)) + EF(j);
               XFORCE2(IDOF(j)) = XFORCE2(IDOF(j)) + EF2(j);
               for k = 1:4
                  XK(IDOF(j), IDOF(k)) = XK(IDOF(j), IDOF(k)) + XKE(j, k);  % Stiffness Assembly
                  XKG(IDOF(j), IDOF(k)) = XKG(IDOF(j), IDOF(k)) + XKGE(j, k);
               end
            end 
        end

        % Apply boundary conditions
        XFORCE(2 * NC - 1) = XFORCE(2 * NC - 1) + FX;
        IBC = [1 2 2*NN];

        if ~isempty(CN)
            IBC = [IBC, CN];
        end

        for i = 1:length(IBC)
             XK(IBC(i), :) = 0;  
             XK(:, IBC(i)) = 0;  
             XK(IBC(i), IBC(i)) = 1; 
             XFORCE(IBC(i)) = 0.0; 
             XFORCE2(IBC(i)) = 0.0;
             XKG(IBC(i), :) = 0;  
             XKG(:, IBC(i)) = 0; 
             XKG(IBC(i), IBC(i)) = 1;
        end

        % Calculate DFORCE and Error
        DFORCE = XFORCE - XFORCE2;
        error = norm(DFORCE);

        % Solve Displacement
        DISP = (XK + XKG) \ (XFORCE - XFORCE2);

           
        k = 1;
        CN = [];

        % Update Positions
        for i = 2:NN
           XG(i) = XG(i) + DISP(2 * i - 1);
           ZG(i) = ZG(i) + DISP(2 * i);
           LENGTH_n(i) = sqrt((XG(i) - XG(i-1))^2 + (ZG(i) - ZG(i-1))^2);
           if ZG(i) <= 0.0 && i>=2
                    ZG(i) = 0.0;
                    XG(i) = XG(i-1) + LENGTH_n(i);
                    CN(k) = 2 * i; % Add element to CN
                    k = k + 1;
           end
        end

        iter = iter + 1;
        
        end
 
        % Plot deformed shape
plot(XG, ZG, 'k', 'LineWidth', 2);
xlabel('Horizontal Position');
ylabel('Vertical Position');
title('Deformation of Cable');
legend('Initial', 'Deformed');
hold on;
end




% % Mirror the deformed cable
% XG_m = 2* XG(NN) -XG; % Reflect the X-coordinates
% ZG_m = ZG;  % Keep the Z-coordinates unchanged
% 
% % Plot the mirrored cable
% plot(XG_m, ZG_m, 'r--', 'LineWidth', 2); % Plot mirrored cable in red dashed line
% 
% % Initial Element Length Calculation
% for i = 1:NE
%     XG_m1 = XG_m(i); XG_m2 = XG_m(i+1);
%     ZG_m1 = ZG_m(i); ZG_m2 = ZG_m(i+1);
%     L0_m(i) = sqrt((XG_m1 - XG_m2)^2 + (ZG_m1 - ZG_m2)^2);
% end




