% LOAD the trace
data=load("RESLT/trace.dat")
% Get the gradient and intercept
m=polyfit(log10(data(:,1)),log10(data(:,2)),1)
% Brewer blue
linestyle = linspecer(1,"qualitative")
% Begin a figure
figure;
hold on
xlabel 'log_{10}h_n'
ylabel 'log_{10} l^2'
% Plot the line and the fit
plot(log10(data(:,1)),log10(data(:,2)),'.','markerfacecolor',linestyle)
plot(log10(data(:,1)),m(1)*log10(data(:,1))+m(2),'-','color',linestyle)
