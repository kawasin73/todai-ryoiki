% load csv file
csv_file = "cputime.csv";
data = csvread(csv_file);

% extract data
% data of dof > 200
extracted = data(9:end,1:4);
dofs = extracted(:,1);
ctimes = extracted(:,2);
ptimes = extracted(:,3);
ebetimes = extracted(:,4);

% convert to log scale
ldofs = log(dofs);
lctimes = log(ctimes);
lptimes = log(ptimes);
lebetimes = log(ebetimes);

% Linear Regression
% URL: https://jp.mathworks.com/help/matlab/data_analysis/linear-regression.html
X = [ones(length(ldofs),1) ldofs];
bc = X\lctimes
be = X\lebetimes

% estimated time in log scale
exlctimes = X * bc;
exlebetimes = X * be;

% estimated time
exctimes = exp(exlctimes);
exebetimes = exp(exlebetimes);

% create graph
scatter(dofs,ctimes);
hold on
scatter(dofs,ebetimes);
loglog(dofs, exctimes);
loglog(dofs, exebetimes);

set(gca,'xscale','log');
set(gca,'yscale','log');
legend('CG method plot', 'PCG with EBE plot', 'CG regression line', 'PCG with EBE regression line', "location", "northwest");
title("calculation time (DOF > 200)");
ylabel("calculation time");
xlabel("degrees of freedom");
grid on
print("cputime-estimate.png", "-dpng")

% estimate the time and DOF
x = (bc(1)-be(1)) / (be(2)-bc(2));
expectedDOF = exp(x)
y = bc(1) + bc(2)*x;
expectedTime = exp(y)
