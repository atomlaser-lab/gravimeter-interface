Irat = [];
delta_0 = [];

figure(6);clf
hold on

% % % Irat 1
x = unique(Irat_1.data.Param);
y = Irat_1.data.R(:,1);
scatter(x,y)
X = 1;
Centre = -165;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % Irat 2
x = unique(Irat_2.data.Param);
y = Irat_2.data.R(:,1);
scatter(x,y)
X = 2;
Centre = -88;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % Irat 4
x = unique(Irat_4.data.Param);
y = Irat_4.data.R(:,1);
scatter(x,y)
X = 4;
Centre = -34;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % Irat 6
x = unique(Irat_6.data.Param);
y = Irat_6.data.R(:,1);
scatter(x,y)
X = 6;
Centre = -7;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % Irat 7
x = unique(Irat_7.data.Param);
y = Irat_7.data.R(:,1);
scatter(x,y)
X = 7;
Centre = -0.015;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % Irat 8
x = unique(Irat_8.data.Param);
y = Irat_8.data.R(:,1);
scatter(x,y)
X = 8;
Centre = 3.5;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% Irat 1/2
x = unique(Irat_1on2.data.Param);
y = Irat_1on2.data.R(:,1);
scatter(x,y)
X = 1/2;
Centre = -236;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % Irat 1/4
x = unique(Irat_1on4.data.Param);
y = Irat_1on4.data.R(:,1);
scatter(x,y)
X = 1/4;
Centre = -300;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % % % Irat 1/6
x = unique(Irat_1on6.data.Param);
y = Irat_1on6.data.R(:,1);
scatter(x,y)
X = 1/6;
Centre = -320;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % % Irat 1/8
x = unique(Irat_1on8.data.Param);
y = Irat_1on8.data.R(:,1);
scatter(x,y)
X = 1/8;
Centre = -330;
Irat = [Irat X];
delta_0 = [delta_0 Centre];

% % % % Irat 1/10
x = unique(Irat_1on10.data.Param);
y = Irat_1on10.data.R(:,1);
scatter(x,y)
X = 1/10;
Centre = -350;
Irat = [Irat X];
delta_0 = [delta_0 Centre];


clear Centre X x y


figure(12);clf
scatter(Irat, delta_0)
xlabel('Intensity Ratio','Interpreter','latex','FontSize',16)
ylabel('Two Photon Detuning (kHz)','Interpreter','latex','FontSize',16)