%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Graphene_5NNFC_3D                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

a = 1.42e-10;
m = 1.99e-26;
pi = 3.141592653589793;
kpoints = 200;

f = diag([25.88, 8.42, 6.183, 4.037, -3.044, -0.492, -3.016, 3.948, 0.516, 0.564, 0.129, -0.521, 1.035, 0.166, 0.110])*16.0217662;

[FA, FB] = rotate(2/3*pi, [1; 0; 0]*a);
[SA, SB] = rotate(1/3*pi, [3/2; sqrt(3)/2; 0]*a);
[TA, TB] = rotate(2/3*pi, [1; sqrt(3); 0]*a);
[LA1, LB1] = rotate(2/3*pi, [2.5; sqrt(3)/2; 0]*a);
[LA2, LB2] = rotate(2/3*pi, [2.5; -sqrt(3)/2; 0]*a);
[WA, WB] = rotate(1/3*pi,[3; 0; 0]*a);
LA = [LA1, LA2];                
LB = [LB1, LB2];

[KAB1, KBA1] = rotate(2/3*pi, f(1:3, 1:3));
[KAA2, KBB2] = rotate(1/3*pi, K(1/6*pi, f(4:6, 4:6)));
[KAB3, KBA3] = rotate(2/3*pi, K(1/3*pi, f(7:9, 7:9)));
[KAB4f, KBA4f] = rotate(2/3*pi, K(acos(2.5/sqrt(7)), f(10:12, 10:12)));
[KAB4s, KBA4s] = rotate(2/3*pi, K(2*pi-acos(2.5/sqrt(7)), f(10:12, 10:12)));
[KAA5, KBB5] = rotate(1/3*pi, f(13:15, 13:15));
KAB4 = cat(3,KAB4f, KAB4s);
KBA4 = cat(3,KBA4f, KBA4s);

kx_path = linspace(-2*pi/(a*3),2*pi/(a*3),kpoints);
ky_path = linspace(-4*pi/(a*3*sqrt(3)),4*pi/(a*3*sqrt(3)),kpoints);

[X,Y] = meshgrid(kx_path,ky_path);

for ii = 1:kpoints
    for ll = 1:kpoints
        if Y(ii,ll) < 1/sqrt(3)*X(ii,ll)-4*pi/(a*3*sqrt(3)) || Y(ii,ll) < -1/sqrt(3)*X(ii,ll)-4*pi/(a*3*sqrt(3)) ...
                ||Y(ii,ll) > -1/sqrt(3)*X(ii,ll)+4*pi/(a*3*sqrt(3)) || Y(ii,ll) > 1/sqrt(3)*X(ii,ll)+4*pi/(a*3*sqrt(3))
            X(ii,ll) = NaN;
            Y(ii,ll) = NaN;
        end
    end
end

k_list = zeros(3,kpoints,kpoints);

for ii = 1:kpoints
    for ll = 1:kpoints
        k_list(:,ii,ll)  = [X(ii,ll), Y(ii,ll), 0];
    end
end

k_num = kpoints*kpoints;
W1 = zeros(kpoints,kpoints,2);
W2 = zeros(kpoints,kpoints,4);
for ii = 1:kpoints
    for jj = 1:kpoints
        k = k_list(:,ii,jj)';
        D = zeros(6, 6);
        DAAs = zeros(3,3);
        DBBs = zeros(3,3);
        DBAs = zeros(3,3);
        DABs = zeros(3,3);
        for ll = 1:3
            DAAs = DAAs + KAA2(:,:,ll) * exp(1i * k * (-SA(:, ll))) + KAA2(:,:,ll+3) * exp(1i * k * (-SA(:, ll+3)))+ ...
               KAA5(:,:,ll) * exp(1i * k * (-WA(:, ll))) + KAA5(:,:,ll+3) * exp(1i * k * (-WA(:, ll+3)));
            DBBs = DBBs + KBB2(:,:,ll) * exp(1i * k * (-SB(:, ll))) + KBB2(:,:,ll+3) * exp(1i * k * (-SB(:, ll+3))) + ...
               KBB5(:,:,ll) * exp(1i * k * (-WB(:, ll))) + KBB5(:,:,ll+3) * exp(1i * k * (-WB(:, ll+3)));
            DABs = DABs + KAB1(:,:,ll) * exp(1i * k * (-FA(:, ll))) + KAB3(:,:,ll) * exp(1i * k * (-TA(:, ll))) +...
                   KAB4(:,:,ll) * exp(1i * k * (-LA(:, ll))) + KAB4(:,:,ll + 3) * exp(1i * k * (-LA(:, ll + 3)));
            DBAs = DBAs + KBA1(:,:,ll) * exp(1i * k * (-FB(:, ll))) + KBA3(:,:,ll) * exp(1i * k * (-TB(:, ll))) + ...
                   KBA4(:,:,ll) * exp(1i * k * (-LB(:, ll))) + KBA4(:,:,ll + 3) * exp(1i * k * (-LB(:, ll + 3)));
        end
        D(1:3, 4:6) = -DABs;
        D(4:6, 1:3) = -DBAs;
        D(1:3, 1:3) = sum(KAB1,3) + sum(KAA2,3) + sum(KAB3,3) + sum(KAB4f,3) + sum(KAB4s,3) + sum(KAA5,3) - DAAs;
        D(4:6, 4:6) = sum(KAB1,3) + sum(KAA2,3) + sum(KAB3,3) + sum(KAB4f,3) + sum(KAB4s,3) + sum(KAA5,3)- DBBs;
        
        D_out = [D(3,3),D(3,6);D(6,3),D(6,6)];
        e1 = eig(D_out);
        w1 = sort(e1);
        W1(ii,jj,:) = real(sqrt(w1))/sqrt(m);

        D_in = [D(1:2,1:2),D(1:2,4:5);D(4:5,1:2),D(4:5,4:5)];
        e2 = eig(D_in);
        w2 = sort(e2);
        W2(ii,jj,:) = real(sqrt(w2))/sqrt(m);
    end
end

W = cat(3,W1(:,:,1),W2(:,:,1),W2(:,:,2),W1(:,:,2),W2(:,:,3),W2(:,:,4))/1e14;

mode = ["ZA";"TA";"LA";"ZO";"TO";"LO"];

d = 4*pi/(3*sqrt(3));
vertices = [0, d; d*sqrt(3)/2, d/2; d*sqrt(3)/2, -d/2;    
            0, -d; -d*sqrt(3)/2, -d/2; -d*sqrt(3)/2, d/2; 0, d];                

figure('OuterPosition',[100 100 1200 800])   
for i=1:6
    subplot(2, 3, i);
    surfc(X*a,Y*a,W(:,:,i),'EdgeAlpha',0.1),shading flat;
    title(mode(i));
    set(gca,'linewidth',1.5,'FontSize',14);
    zlim([0 3.1]);
    xlim([-2*pi/3-0.2,2*pi/3+0.2]);
    ylim([-4*pi/(3*sqrt(3))-0.2,4*pi/(3*sqrt(3))+0.2]);
    xlabel('$k_x a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    ylabel('$k_y a$','Interpreter','latex','FontSize',20,'FontWeight','bold');
    zlabel('$ \omega,10^{14} $ rad/s','Interpreter','latex','FontSize',20,'FontWeight','bold');
    hold on
    plot(vertices(:, 1), vertices(:, 2), 'b-');
end

save data.mat X Y W

function [R, RB] = rotate(theta, r)
    Um = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    inv_Um = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    U = [cos(pi), sin(pi), 0; -sin(pi), cos(pi), 0; 0, 0, 1];
    inv_U = [cos(pi), -sin(pi), 0; sin(pi), cos(pi), 0; 0, 0, 1];
    n = 2*pi/theta;
    if numel(r) == 3
        R = zeros(3,n);      
        RB = zeros(3,n);
    else
        R = zeros(3,3,n);      
        RB = zeros(3,3,n);
    end
    for i = 1:n
        if numel(r) == 3
            r = inv_Um * r;
            rb = r*(-1);
            R(:, i) = r;
            RB(:, i) = rb;
        else
            r = inv_Um * r * Um;
            rb = inv_U * r * U;
            R(:, :, i) = r;
            RB(:, :, i) = rb;
        end

    end
end

function k = K(theta, k)
    U = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    inv_U = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    k = inv_U * k * U;
end

