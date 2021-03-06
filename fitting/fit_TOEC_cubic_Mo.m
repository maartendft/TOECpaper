clear all; clc; close all; format long;
load D1_Mo; load D2_Mo; load D3_Mo; load D4_Mo; load D5_Mo; load D6_Mo

s1 = 2;

%%
% deformation mode 1, t_1
figure(1)
hold on
plot(D1_Mo(:,1), D1_Mo(:,2), '.')
f1 = polyfit(D1_Mo(s1:end,1), D1_Mo(s1:end,2)/10, 2);
plot(D1_Mo(:,1), 10*polyval(f1, D1_Mo(:,1)))

C11 = f1(2);
C111 = f1(1)*2;
C11
C111


%%
% deformation mode 1, t_2
figure(2)
plot(D1_Mo(:,1), D1_Mo(:,3)/10, '.')
f2 = polyfit(D1_Mo(s1:end,1), D1_Mo(s1:end,3)/10, 2);
C12 = f2(2);
C112 = f2(1)*2;
C12
C112


%%
% deformation mode 2, t_1
figure(3)
plot(D2_Mo(:,1), D2_Mo(:,2)/10, '.')
f3 = polyfit(D2_Mo(s1:end,1), D2_Mo(s1:end,2)/10, 2);
f3
% C12 = f3(2)/2
C123 = f3(1)-C112


%%
% deformation mode 3, t_4
figure(4)
plot(D3_Mo(:,1), D3_Mo(:,2)/10, '.')
f4 = polyfit(D3_Mo(6:15,1), D3_Mo(6:15,2)/10, 1);
f4;
C44 = f4(1)/2


%%
% deformation mode 4, t_4
figure(5)
plot(D4_Mo(:,1), D4_Mo(:,2)/10, '.')

f5 = polyfit(D4_Mo(s1:end,1), D4_Mo(s1:end,2)/10, 2);
f5

% C44 = f5(2)/2
C144 = f5(1)/2

%%
% deformation mode 5, t_5
figure(6)
plot(D5_Mo(:,1), D5_Mo(:,2)/10, '.')

f6 = polyfit(D5_Mo(s1:end-3,1), D5_Mo(s1:end-3,2)/10, 2);
f6

% C44 = f6(2)/2
C155 = f6(1)/2


%%
% deformation mode 6, t_4
figure(7)
plot(D6_Mo(:,1), D6_Mo(:,2)/10, '.')

f7 = polyfit(D6_Mo(6:12,1), D6_Mo(6:12,2)/10, 2);
f7

C44 = f7(2)/2
C456 = f7(1)/4




