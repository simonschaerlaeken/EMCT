% Quian Quiroga's method for detecting leadership
% this script implements the Quian Quiroga's method to measure Event
% Synchronicity (ES)and time delay patterns (Phys.Rev. E 66 041904). It can
% be apply to univariate time series or
% multivariate time series rearranged in probability domain.


%input: peakfinal, t x 2 matrix containing in each row time and peak value
%       t, length of the input vectors
%       delta n

%output: if 'tot' is specified, the ouptut is the pair Qtau and qtau
%        measuring the synchonization and the delay behavior between the
%        whole signals;
%        if 'tsl' is specified, the output are the time resolved variants
%        of Qtau and qtau;
%        if 'ssl' is specified, the output are Qtau and qtau averaged
%        over the last delta_n time steps.



function [Qout,qout]= Qqm2(varargin)

error(nargchk(1,8,nargin));
if nargout>2, error('Too many output arguments'), end

check_type={'tot','tsl','ssl'};	% total, time resolved, sample resolved


[s1,s2]=size(varargin{1});
[s3,s4]=size(varargin{2});

    if ((s2~=2)|(s4~=2))
       error('too many columns in the input vectors! the columns must be 2');
    else
       peakfinal1=varargin{1};
       peakfinal2=varargin{2};
    end
    
%checking type
i_char=find(cellfun('isclass',varargin,'char'));
temp_type=0;
  if ~isempty(i_char)
    for i=1:length(i_char), 
         varargin{i_char(i)}(4)='0';
         temp_type=temp_type+strcmpi(varargin{i_char(i)}(1:3),check_type'); 
    end
    type=min(find(temp_type));
      if isempty(type), type=2; end
      if type>3, type=3; end
  end

if type==3
    delta_n=varargin{5};
end

make_plot = 0;
if nargin == 8
    make_plot = varargin{8};
end
    
t=max(max(peakfinal1(:,1)),max(peakfinal2(:,1)));
if nargin >= 7
    dim = varargin{7};
    if dim > t
        t = dim;
    end
end

% pulses-train approximation of the input time series called peakfinal 
ts1=zeros(1,t);
peakfinal1=peakfinal1(:,1);
a=find(peakfinal1==0);
peakfinal1(a)=[];

for num=1:1:size(peakfinal1,1)
    ts1(1,peakfinal1(num))=1;
end

ts2=zeros(1,t);
peakfinal2=peakfinal2(:,1);
b=find(peakfinal2==0);
peakfinal2(b)=[];


for num=1:1:size(peakfinal2,1)
    ts2(1,peakfinal2(num))=1;
end



% choice of optimal time lag tau
delta1(1)=peakfinal1(1);
for count1=1:1:(size(peakfinal1,1)-1)
    delta1(count1+1)= peakfinal1(count1+1)-peakfinal1(count1);
end

delta2(1)=peakfinal2(1);
for count2=1:1:(size(peakfinal2,1)-1)
   delta2(count2+1)= peakfinal2(count2+1)-peakfinal2(count2);
end
   delta_min=min(min(delta1),min(delta2))/2;
   %delta_min=3;



switch(type)
  case{1} %computation of Q(tau) and q(tau)
%c(x/y)tau
for num_peak_x=1:1:size(peakfinal1,1)
    jay_out_y=0;
    for num_peak_y=1:1:size(peakfinal2,1)
        diff=(peakfinal1(num_peak_x)-peakfinal2(num_peak_y));
        if ((diff>0) & (diff<=delta_min))
            jay_out=1;
        elseif diff==0
            jay_out=1/2;
        else
            jay_out=0;
        end
          jay_out_y=jay_out_y+jay_out;
          
    end
    jay_out_xy(num_peak_x)=jay_out_y;
end
cxy=sum(jay_out_xy);

%c(y/x)tau
for num_peak_y=1:1:size(peakfinal2,1)
    jay_out_x=0;
    for num_peak_x=1:1:size(peakfinal1,1)
        diff=(peakfinal2(num_peak_y)-peakfinal1(num_peak_x));
        if ((diff>0) & (diff<=delta_min))
            jay_out=1;
        elseif diff==0
            jay_out=1/2;
        else
            jay_out=0;
        end
          jay_out_x=jay_out_x+jay_out;
          
    end
    jay_out_yx(num_peak_y)=jay_out_x;
end

 cyx=sum(jay_out_yx);

Q_tau= ((cxy + cyx)/sqrt(size(peakfinal1,1)*size(peakfinal2,1)));
q_tau= ((cyx - cxy)/sqrt(size(peakfinal1,1)*size(peakfinal2,1)));

Qout=Q_tau;
qout=q_tau;




case{2} %computation of Q(n) and q(n): time-resolved version of Q(tau) and q(tau)
%cn(x/y)
jay_fincnxy=zeros(1,t);
for num_peak_x=1:1:size(peakfinal1,1)
    jay_out_y=0;
    for num_peak_y=1:1:size(peakfinal2,1)
        diff=(peakfinal1(num_peak_x)-peakfinal2(num_peak_y));
        if ((diff>0) & (diff<=delta_min))
            jay_out=1;
        elseif diff==0
            jay_out=1/2;
        else
            jay_out=0;
        end
          jay_out_y=jay_out_y+jay_out;

          for sample=1:1:t
              jay(sample)=jay_out_y*heaviside(sample-peakfinal1(num_peak_x));
              if(isnan(jay(sample)))
                  jay(sample)=0;
              end
          end
    end
    jay_fincnxy=jay_fincnxy+jay;
end

%cn(y/x)
jay_fincnyx=zeros(1,t);
for num_peak_y=1:1:size(peakfinal2,1)
    jay_out_x=0;
    for num_peak_x=1:1:size(peakfinal1,1)
        diff=(peakfinal2(num_peak_y)-peakfinal1(num_peak_x));
        if ((diff>0) & (diff<=delta_min))
            jay_out=1;
        elseif diff==0
            jay_out=1/2;
        else
            jay_out=0;
        end
          jay_out_x=jay_out_x+jay_out;

          for sample=1:1:t
              jay(sample)=jay_out_x*heaviside(sample-peakfinal2(num_peak_y));
              if(isnan(jay(sample)))
                  jay(sample)=0;
              end
          end
    end
    jay_fincnyx=jay_fincnyx+jay;
end


Qn=jay_fincnxy+jay_fincnyx;
qn=jay_fincnyx-jay_fincnxy;

Qout=Qn;
qout=qn;

subplot(1,2,1);
plot(Qout)
xlabel('samples');
title('Qn');

subplot(1,2,2);
plot(qout)
xlabel('samples');
title('qn');

case{3} %computation of Q'(n) and q'(n): averaged variants of Q(n) and q(n) over the last delta_n time steps
     
if ((delta_n<delta1) & (delta_n<delta2))
   delta_n=delta1+delta2;
   disp('too small deltaT, the algorithm is using a default value');
end

%jx
jay_x = zeros(size(peakfinal1, 1), 1);
for num_peak_x = 1 : 1 : size(peakfinal1, 1)
    for num_peak_y = 1 : 1 : size(peakfinal2, 1)
        diff = peakfinal1(num_peak_x) - peakfinal2(num_peak_y);
        if ((diff > 0) && (diff <= delta_min))
            curr_jay = 1;
        elseif diff == 0
            curr_jay = 1/2;
        else
            curr_jay= 0;
        end
        jay_x(num_peak_x) = jay_x(num_peak_x) + curr_jay;
    end
end

%jy
jay_y = zeros(size(peakfinal2, 1), 1);
for num_peak_y = 1 : 1 : size(peakfinal2, 1)
    for num_peak_x = 1 : 1 : size(peakfinal1, 1)
        diff = peakfinal2(num_peak_y) - peakfinal1(num_peak_x);
        if ((diff > 0) && (diff <= delta_min))
            curr_jay = 1;
        elseif diff == 0
            curr_jay = 1/2;
        else
            curr_jay = 0;
        end
          jay_y(num_peak_y) = jay_y(num_peak_y) + curr_jay;
    end
end

Qpn = zeros(1, t);
qpn = zeros(1, t);
Qpn(1 : delta_n) = nan;
qpn(1 : delta_n) = nan;

for n = delta_n + 1 : 1 : t
    theta_n_x = heaviside(n - peakfinal1);
    theta_n_y = heaviside(n - peakfinal2);
    theta_dn_x = heaviside(n - delta_n - peakfinal1);
    theta_dn_y = heaviside(n - delta_n - peakfinal2);
    deltanx = size(find(ts1((n - delta_n) : n)),2); %number of events in [n-deltan,n] for x
    deltany = size(find(ts2((n - delta_n) : n)),2); %number of events in [n-deltan,n] for y]
    if (deltanx == 0 || deltany == 0)
        if ((theta_n_y - theta_dn_y)' * jay_y + (theta_n_x - theta_dn_x)' * jay_x) == 0
            Qpn(n) = 0;
        else
            Qpn(n) = nan;
        end 
        if ((theta_n_y - theta_dn_y)' * jay_y - (theta_n_x - theta_dn_x)' * jay_x) == 0
            qpn(n) = 0;
        else
            qpn(n) = nan;
        end
    else
        Qpn(n) = ((theta_n_y - theta_dn_y)' * jay_y + (theta_n_x - theta_dn_x)' * jay_x) / sqrt(deltanx * deltany);
        qpn(n) = ((theta_n_y - theta_dn_y)' * jay_y - (theta_n_x - theta_dn_x)' * jay_x) / sqrt(deltanx * deltany); 
    end
end

Qout=Qpn;
qout=qpn;

if make_plot
    subplot(1,2,1);
    plot(Qout)
    xlabel('samples');
    title('Averaged Qn');

    subplot(1,2,2);
    plot(qout)
    xlabel('samples');
    title('Averaged qn');
end

end







