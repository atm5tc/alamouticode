% Alamouti Coding & Demodulation
clear all;
load 'alamouti.mat';
h11 = 0.8*exp(j*pi/8);
h12 = 0.5*exp(-j*pi*3/4);
%% Create Scatter Plot of 'data' and 'r'
figure(1);
hold on;
plot(real(data),imag(data),'ob');
plot(real(r),imag(r),'dr');
legend('data','r');
xlabel('Real Part Re[point]');

ylabel('Imag Part Im[point]');
title('Scatter Plot of Data & Received Signal');
hold off;
%close all;
%% Demodulate & Decide on symbols
%  Convert to 2x2 matrix, in reference to 2 channel demodulation scheme
close all;
% Turn data into 2x(data/2) matrix
x = [r(1:2:50);r(2:2:50)];
y = zeros(size(x));
%% Demodulate rx with Alamouti Code
for i=1:length(x)
    y(1,i) = conj(h11)*x(1,i)+h12*conj(x(2,i));
    y(2,i) = 1*conj(h12)*x(1,i)-1*h11*(conj(x(2,i)));
end
%% Process symbols according using square method
z = zeros(1,size(y,1)*size(y,2));
z(1:2:end) = (abs(h11)^2+abs(h12)^2)*y(1:2:end);
z(2:2:end) = (abs(h11)^2+abs(h12)^2)*y(2:2:end);
%% Complex plot of reception post demod
hold on;
plot(real(z),imag(z),'+');
xlabel('Real Part Re[point]');
ylabel('Imag Part Im[point]');
title('Scatter Plot of Demodulated Symbols');
hold off;
%% Determine phase rotation by finding closest point to 0
theta = exp(j*pi);
for i=1:length(z)
    if(real(z(i))>=0 && imag(z(i)) >=0)
        if(angle(z(i))<angle(theta))
            theta = z(i);
        end
    end
end
theta = angle(theta);
for i=1:length(z)
    z(i) = mag(z(i))*exp(angle(z(i))*j+theta*j);
end
%% Make Decisions on Calculated Data Received
dec = 8; % 2^n # of decisions to make, for QPSK this is 4, for us it's 8
r_data = zeros(dec,size(y,1)*size(y,2)); % Matrix to calculate the 'closeness' of each 'dec' decision
for i=1:length(z)
    for d = 1:dec
        r_data(d,i) = abs((real(z(i))-real(exp((d-1)*j*pi/(dec/2)))))+...
             j*abs(imag(z(i))-imag(exp((d-1)*j*pi/(dec/2))));
    end
end
%
%
r = [r;zeros(size(z))];
[ans,r(2,:)] = min(abs(r_data));
r(2,:) = r(2,:)-1;
% Note at this point rx'd demodulated symbols are stored in r(col 1)
%% Decisions on Original Data Set
dec = 8;
data = [data;zeros(dec,length(data))];
%
for k=1:length(data)
    for d = 1:dec
        data(d+1,k) = abs((real(data(1,k))-real(exp((d-1)*j*pi/(dec/2)))))+...
             j*abs(imag(data(1,k))-imag(exp((d-1)*j*pi/(dec/2))));
    end
end
[ans,r(1,:)] = min(data(2:dec+1,:));
r(1,:) = r(1,:)-1;

%% Create stem plot of data correlation
stem(r(1,:)-r(2,:));
xlabel('Symbol Index');
ylabel('Correlation [data_{rx}-data_{tx}]');
title('Stem Plot of Demodulated Symbols');
for i=1:length(x)
    y(1,i) = h11*x(1,i)+h12*x(2,i);
    y(2,i) = h11*(-1*conj(x(2,i)))+h12*conj(x(1,i));
end



%% Make decisions on received array
dec = 8; % 2^n # of decisions to make, for QPSK this is 4, for us it's 8
z = zeros(dec,size(y,1)*size(y,2)); % Matrix to calculate the 'closeness' of each dec decisions
data_r = zeros(1,size(y,1)*size(y,2));

% Decisions on Calculated Data Received
for i=1:length(z)
    for d = 1:dec
        z(d,i) = abs((real(y(i))-real(exp((d-1)*j*pi/(dec/2)))))+...
             j*abs(imag(y(i))-imag(exp((d-1)*j*pi/(dec/2))));
    end
    [data_r(1,i),data_r(2,i)] = min(z(:,i));
end
data_r(2,:) = data_r(2,:)-1;
%%
% Decisions on Original Data Set
dec = 8;
data = [data;zeros(dec,length(data))];

for k=1:length(data)
    for d = 1:dec
        data(d+1,k) = abs((real(data(k))-real(exp((d-1)*j*pi/(dec/2)))))+...
             j*abs(imag(data(k))-imag(exp((d-1)*j*pi/(dec/2))));
    end
    [ans,data(1,k)] = min(data(2:dec+1,k));
end


%%
%plot(real(z),imag(z),'ob');
stem(abs(data_r(2,:)-data_r(4,:)));