%%
% SECTION A
t = 0:6;
t = 60*60*24 .* t; t= t.'

v = [12001.9 12010.9 12018.1 12027.2 12036.9 12044.1 12048.0]

xhi = @(k,j) sum(t.^(k+j));

Xhi = [xhi(0,0) xhi(0,1);
       xhi(1,0) xhi(1,1)]
gamma = @(k) sum(t.^k .* v.');
Gamma = [gamma(0); 
         gamma(1)]

a = Xhi^(-1) * Gamma
ac = a(2); v0 = a(1);
fprintf('V(t) =  %d t + %d [m/s]\n', [ac v0])


%%
% SECTION B
F = ac * 500


%%
% SECTION C

fv = ac*t + v0
y_bar = mean(v)
ss_tot = sum((v-y_bar).^2)
ss_res = sum((fv-y_bar).^2)
Rsqure = 1 - ss_res/ss_tot

%%
% FOR MYSELF
% verify by visualation the results of the Q above..
% 
% T = linspace(0,t(end))
% plot(T, ac*T+v0)
% title 'visual verify...'
% ylabel 'V(t) [m/s]'
% xlabel 't[sec]'
% hold on
% plot(t, v, 'rX', MarkerSize=10,LineWidth=3)
% hold off
