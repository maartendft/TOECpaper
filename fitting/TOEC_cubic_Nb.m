clear all; clc; close all; format long;
load D1_Nb; load D2_Nb; load D3_Nb; load D4_Nb; load D5_Nb; load D6_Nb
%%
% deformation mode 1, t_1
s1=4;
figure(1)
hold on
plot(D1_Nb(4:end,1), D1_Nb(4:end,2)/10, '.')
f1 = polyfit(D1_Nb(s1:end,1), D1_Nb(s1:end,2)/10, 2);
plot(D1_Nb(4:end,1), 1*polyval(f1, D1_Nb(4:end,1)))
xlabel('\eta_1')
ylabel('t_1')
C11 = f1(2);
C111 = f1(1)*2;
C11
C111
print -dpng -r400 pic1

%%
% deformation mode 1, t_2
s1=2;
figure(2)
hold on
plot(D1_Nb(2:end,1), D1_Nb(2:end,3)/10, '.')
f2 = polyfit(D1_Nb(s1:end,1), D1_Nb(s1:end,3)/10, 2);
plot(D1_Nb(2:end,1), polyval(f2, D1_Nb(2:end,1)))
xlabel('\eta_1')
ylabel('t_2')
C12 = f2(2);
C112 = f2(1)*2;
C12
C112
print -dpng -r400 pic2

%%
% deformation mode 2, t_1
figure(3)
hold on
plot(D2_Nb(2:end,1), D2_Nb(2:end,2)/10, '.')
f3 = polyfit(D2_Nb(s1:end,1), D2_Nb(s1:end,2)/10, 2);
plot(D1_Nb(2:end,1), polyval(f3, D1_Nb(2:end,1)))
xlabel('\eta_2 = \eta_3')
ylabel('t_1')

f3
% C12 = f3(2)/2
C123 = f3(1)-C112
print -dpng -r400 pic3

%%
% deformation mode 3, t_4
figure(4)
plot(D3_Nb(:,1), D3_Nb(:,2)/10, '.')
f4 = polyfit(D3_Nb(6:12,1), D3_Nb(6:12,2)/10, 1);
xlabel('\eta_4')
ylabel('t_4')
f4;
C44 = f4(1)/2
print -dpng -r400 pic4

%%
% deformation mode 4, t_4
figure(5)
hold on
plot(D4_Nb(:,1), D4_Nb(:,2)/10, '.')
f5 = polyfit(D4_Nb(4:end,1), D4_Nb(4:end,2)/10, 2);
plot(D1_Nb(:,1), polyval(f5, D1_Nb(:,1)))
xlabel('\eta_1 = \eta_4')
ylabel('t_4')
f5

f5(2)/2
C144 = f5(1)/2
print -dpng -r400 pic5
%%
% deformation mode 5, t_5
figure(6)
hold on
plot(D5_Nb(:,1), D5_Nb(:,2)/10, '.')
f6 = polyfit(D5_Nb(5:end,1), D5_Nb(5:end,2)/10, 2);
plot(D1_Nb(:,1), polyval(f6, D1_Nb(:,1)))
xlabel('\eta_1 = \eta_5')
ylabel('t_5')

f6

f6(2)/2
C155 = f6(1)/2
print -dpng -r400 pic6

%%
% deformation mode 6, t_4
figure(7)
hold on
plot(D6_Nb(:,1), D6_Nb(:,2)/10, '.')
f7 = polyfit(D6_Nb(7:end-5,1), D6_Nb(7:end-5,2)/10, 2);
plot(D1_Nb(:,1), polyval(f7, D1_Nb(:,1)))
xlabel('\eta_4 = \eta_5 = \eta_6')
ylabel('t_4')

f7

f7(2)/2
C456 = f7(1)/4
print -dpng -r400 pic7



