nntwarn off;

p1 = [1 1000 0 15;1 800 200 15;1 600 400 15;1 400 600 15;...
3 1000 0 15;3 800 200 15;3 600 400 15;3 400 600 15;...
10 1000 0 15;10 800 200 15;10 600 400 15;10 400 600 15;...
30 1000 0 15;30 800 200 15;30 600 400 15;30 400 600 15;...
100 1000 0 15;100 800 200 15;100 600 400 15;100 400 600 15;...
365 1000 0 15;365 800 200 15;365 600 400 15;365 400 600 15]';

p2 = [1 1000 0 100;1 800 200 100;1 600 400 100;1 400 600 100;...
3 1000 0 100;3 800 200 100;3 600 400 100;3 400 600 100;...
10 1000 0 100;10 800 200 100;10 600 400 100;10 400 600 100;...
30 1000 0 100;30 800 200 100;30 600 400 100;30 400 600 100;...
100 1000 0 100;100 800 200 100;100 600 400 100;100 400 600 100;...
365 1000 0 100;365 800 200 100;365 600 400 100;365 400 600 100]';

p = [p1, p2];

t1 = [0.119;1.5303;3.105;4.4253;...
0.5113;0.5993;1.697;3.0333;...
0.675;0.8683;1.004;1.0973;...
0.2737;0.183;0.2287;0.3017;...
0.2593;0.224;0.1973;0.208;...
0.011;0.0103;0.013;0.023]';

t2 = [5.3;68;138;196.7;...
22.7333;26.6;75.4;134.8;...
30;38.5667;44.6;48.7333;...
12.2;8.1333;10.1667;13.3667;...
11.5333;9.9667;8.8;9.2333;...
0.5333;0.5;0.5667;1.0333]';

t = [t1, t2];

[pn,ps] = mapminmax(p);%normalizing input data
[tn,ts] = mapminmax(t);%normalizing output data
net=newff(minmax(pn),[9 1],{'logsig','purelin'},'trainlm');%creat bp neural network
% net=feedforwardnet(7);%creat neural network and the hidden layer number is 5
% set parameters 
% net.layers{1}.transferFcn='tansig';%from input layer to hidden layer transfer function
% net.layers{2}.transferFcn='purelin';%from hidden layer to output layer transfer function
%net.trainFcn='traingdx';%Gradient descent with momentum and adaptive learning rate backpropagation
% net.trainFcn='trainlm';
% net.performFcn='mse';
net.trainParam.epochs=10000;
net.trainParam.mc=0.9;%Momentum constant
net.trainParam.lr=0.01;%Learning rate
net.trainParam.goal= 0.0001;
net=init(net);%Initializes the network
net1=train(net,pn,tn);%train the neural network
an=sim(net1,pn);%simulate bp neural network
a = mapminmax('reverse',an,ts);
performance = perform(net1,t,a);%mean squared error
RMSE = sqrt(performance);%root mean squared error

figure(1)
plot(a,':og');  
hold on  
plot(t,'-*');  
legend('Predicted output','Actual output');



