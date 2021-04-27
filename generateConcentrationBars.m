close all
clc
clear all

Pset=[0.0196161941252548,0.0845081899260819,0.00854050392445480,23.9720980151291,0.148031000203075,0.634874372630089,0.128402833632652,0,20.0540479883642,3.19597545909701,0.172667672916842,0.0348848073266119,0.000400366818540836,0.152618607413662,0.607898902070148,0.0157090661451648,0.626784774174197,0.0925032952587944,0.0293013990558697,0.203202811615799,0.229731532848886,0.000137666970461378,0.00203215228242153,0.000442421489411307,0.841140753556429,0.153778252055787,0.125526764113365,0.148950619425383,49.7988347631407,237.829687748671,6.16081895526424,0.141158792243207,0.0676568351264791,1.24802676356656,0.334982163601985,0.109687694753112,0.0245882508301444,0,0.390258640480244,2.64946031309526e-06,0.0571220074248525,0.0239123613880403,0.835481986431901,0.546601182386169,0.0397035905821537,0.0100793853545540,0.104704214385541,0.0471279183840516,6.44283937130162,0.00375955556775643,0.00302485389871826,0.256499287452978,0.0298548125643264,12689.0255780500,0.115306309712000,24.9868269193688,0.00212683479605392,0,0.000506859495373619,33.8593998301064,0.000809568509023799,0.825708493125192,0.0130760062278384,0.0256948153839262,0.509838481449771,0.00110186311930334,0.192415146159320,0.112846265577634,3.60453068978040,0.127826376008879,0.0115085419530951,6.53460085845081,0.00344334056239891,0.0673172236089235,0.00724548527098653,0.205711471361743,50.3756851350204,0.0199275940236153,0.182431306193620,0.342063054281516,0.463891508530205,0.0134106031314399,0.00577268116908457,0.0612702365031117,29.8950044749535,0.658693606789406,0.0556373948903302,1.19350218026503,0.113190170599104,0.0203610031928475,0.0776271722560993,25.6263976671721,62.4341831170576,0.00262688690580012,0.0182890755107398,6.46505475724398,0.262809546780673,0.00388975172376345,0.234199350377072,1402.67572918751,65.8034275104182,2.48838852886230,0.0538980291597912,0.0895734337651592,0.00318965765203073,0.0935708379122437,1.58817641734158,0.431255021446841,0.891694257015847,0.0569081596554034,0.0545223860722562,0.0241780728488369,0.00427882790595702,0.156146185701553,0.0252202596567516,0.0314944841881572,1.03607330764189,0.0190055062472176,3.43369064092679,0.158733612598239,0.178255754895477,3.74468599165674,0.162065238187832,0.210175428062283,0];
Pset=[0.0300258560322406,0.0700450577585854,0.0134855256655188,25.6214008379448,0.138883706102963,0.858105905924659,0.0876006037055290,0.0658797679919270,11.2708649987539,3.84253971168942,0.295536837356773,0.0255831882697163,0.000439067681752294,0.201857276256458,0.782440137310070,0.0158221757906231,0.477762481247816,0.0872534143021344,0.0226115875552263,0.357065236163759,0.221860810142983,0.000457053247783196,0.00155432515815765,0.000349879370414712,1.60671551277495,0.109227676918023,0.0877988044664160,0.112738079744207,28.8794785511875,104.397200819046,9.31642514221566,0.140050278945431,0.0619751931560969,1.24802676356656,0.356692896893517,0.0427652762289979,0.0325741748775381,0,0.173257818790992,7.78595288748276e-05,0.0576379224901726,0.0452606206389411,0.181581959664141,0.892245143712699,0.0309463574909633,0.0113746158943527,0.134262880682973,0.0425019683576624,13.8462140729987,0.00124657933912000,0.00236054239713399,0.279678469214574,0.0197127816226931,3572.07812987308,0.0821987559327118,50.5419701826314,0.00841098579954952,0,0.000667711735938746,23.3845894453517,0.000927565158519696,1.01982421944280,0.00744127770137666,0.0319535883877096,0.881115026633867,0.000698843389954051,0.197987904902603,0.0823545657129377,1.76805488255182,0.136900151075001,0.00852844706291415,4.63604054600568,0.00281753089173504,0.0604836986591620,0.00711924382189207,0.204899309574355,39.3379901067497,0.0277299713582054,0.101814854068068,0.625242871692251,0.401461802451721,0.0197854850087685,0.00441196402667265,0.0378391789861344,17.0767855822946,0.375622331721223,0.270134082256646,2.17744417133131,0.0987613680961496,0.0188046600972326,0.0589734930219092,19.8024778752044,80.1855616471112,0.00300500456726273,0.0187772364969706,7.26536435640985,0.704987530864179,0.00248907281722838,0.199012100991402,1165.53508175719,28.5526756493611,1.62060776507696,0.0543554945701113,0.130307251742549,0.0191918573916840,0.0935708379122437,1.43476866642204,0.518583945644741,0.959843057292479,0.0548611179208426,0.0766513023363724,0.122397175316066,0.0155062075983528,0.357232775997917,0.0240071927029249,0.0639083333042296,0.836641564007736,0.0118314121599819,3.17409160243149,0.111226681037340,0.199225999054764,2.37391666848024,0.142483857120125,0.206820072842497,0];
celltype = 'SW';
[tout, teout, ieout, yout, cellCycleIniTimes] = CauloSim(celltype,'WT', 1800,Pset);
graphCellCycle(tout, teout, ieout, yout, cellCycleIniTimes, Pset, 'ST')
[t, y, yayornay, events]=getLastCycle2(tout,yout,teout,ieout);
% 
% ConcentrationBars(t,y,[2,3,4,8],{'CtrA~P','GcrA','DnaA','CcrM'},[1 0 0; 0 1 0; 0 0 1; 1 0 1])

ConcentrationRings(t,y,[2,4,3,8],{'CtrA~P','GcrA','DnaA','CcrM'},[1 0 0; 0 1 0; 0 0 1; 1 0 1], 3.5, 3.7,events)



% 
% 
% colors=[1 0 0];
% xf=0; yf=0;
% r=2;
% R=2.25;
% figure()
% axis equal
% % set(gcf, 'Position', [100 300 1200 300]);
% hold on
% % t= [1 2 3 4 5]
% p = linspace(0,2*pi,length(t));
% j=1;
% scale=1;
% plotvalues=2;
% yset=y;
% for i=1:length(t)-1
%     colorshade=[(1+(colors(j,1)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,2)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale),(1+(colors(j,3)-1)*yset(i,plotvalues(j))/max(yset(:,plotvalues(j))))^(1/scale)];
%     colorshade(colorshade > 1)=1;
%     %     h=fill([p(i), p(i+1), t(i+1), t(i)], ...
%     %         [spacings(j), spacings(j), spacings(j)-1, spacings(j)-1],colorshade);
%     x = r*cos(p(i:i+1));
%     y = r*sin(p(i:i+1));
%     X = R*cos(p(i:i+1));
%     Y = R*sin(p(i:i+1));
%     h = patch([x flip(X)],[y flip(Y)],colorshade);
%     set(h,'EdgeColor','none')
% end
% p = linspace(0,2*pi,length(t)-1);
% x = xf + r*cos(p);
% y = yf + r*sin(p);
% X = xf + R*cos(p);
% Y = yf + R*sin(p);
% % P = patch([x X],[y Y],[1,.5,.5],'linestyle','non','facealph',.5);
% L(1) = line(x,y,'color','k');
% L(2) = line(X,Y,'color','k');