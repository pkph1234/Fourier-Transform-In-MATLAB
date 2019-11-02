clc;
close all;
x=menu('Select the signal to plot','1) Step Function','2) Ramp Function','3) Sinosudial Function','4) Exponential Function');
switch x
case 1
   clc;
   close all;
   disp('This is Step function:')
%% Intilization for making the function:
% min='what is your min range?';
% Min=input(min);
Min=-10;% Intialize the lower boundry for plotting
% max='What is your max range?';
% Max=input(max);
Max=10;% initialize of higher boundry for plotting
% sample='how many sample you want?';
% Sample=input(sample);
Sample=10000;% predefind number of samples
Stepwidth=input('How much is the step width?');
% Stepwidth=1;
a=Sample/(Max-Min);%for sampling terms to find sapcing between each sample
% delay='How much time shift do you need?';
Delay=input('What is your delay in second?');
% Delay=0;
def=Delay-Min;% to shift the axis accoring to delay
Def=Min-Min;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
Amp=input('What is the Amplitude?');
% Amp=1;
s=1/a;
d=0;% Start value number to find the value of A0
T0=Max-Min;
L=T0/2;
As=length(n);
ft=zeros(1,As);
%% making array of F(t) , A0 , An & Bn
for x=1:1:length(n)% making array of F(t)
    if x>=(def)*a && x<=(def+Stepwidth)*a
        y(x)=Amp;
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)% plotting the sin function for the calculation of An & Bn
        z(x)=sin(d*pi*x/(a*L));
    else
        z(x)=0;
        end
        if x>=Def*a && x<=def*a+(Max*a)
        w(x)=cos(d*pi*x/(a*L));
    else
        w(x)=0;
        end
end
%% Value of A0,An,Bn;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn

for x=1:1:length(n)% Plotting the A0 at t=0
if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=A0;   
else
    ft(x)=0;
end
end
%% Starting the loop for adding the terms
for D=1:1:200000% Creating loop for adding the value of n
    
for x=1:1:length(n)   
    if x>=(def)*a && x<=(def+Stepwidth)*a
        y(x)=Amp;
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(D*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(D*pi*x/(a*L));
    else
        w(x)=0;
        end
end
%% Value of A0, An & Bn
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
%% Finding the array of summing term
for x=1:1:length(n)% loop for getting added function each time before n value is incremented by 1

if x>=(Def)*a && x<=(def)*a+(Max*a)
      ft(x)=ft(x)+(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ft(x)=0;
end
end
for x=1:1:length(n)% loop for finding the individual added function at the value of n

if x>=(Def)*a && x<(def)*a+(Max*a)
      ftt(x)=(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ft(x)=0;
end
end
%% Finding error
RMSE=(sqrt(mean((y-ft).^2)))*100;
% asd=(wee-we)*100
if RMSE<=5
    break
end
%% Plotting the graphs
  subplot(3,1,1)
  plot(n,y,'r','linewidth',2) 
  axis auto
  grid on
  xlabel('Time(in second)')
  ylabel('Amplitude')
  title('Input Step Function F(t)')
subplot(3,1,2)
 n=Min:(Max-Min)/Sample:Max;
plot(n,ftt,'linewidth',2)
axis auto
grid on
xlabel('Time (in second)');
ylabel('Amplitude');
  title(['Summing Signal at n=',num2str(D), ', A0=',num2str(A0),', A',num2str(D),'=',num2str(An),', B',num2str(D),'=',num2str(Bn)]);
subplot(3,1,3)
plot(n,ft,'b','linewidth',2)
axis auto
grid on
xlabel('Time(in second)')
  ylabel('Amplitude')
  title(['Output Step Function SF(t)',', Error in the plot is: ',num2str(RMSE),'%'])
pause(1)

end


case 2
     disp('This is ramp Function')  
clc
close all
Min=-10;% intialize the lower boundry
Max=10;% intialize the upper boundry
Sample=10000;% initialize the sample size
a=(Max-Min)/Sample;
End=input('From where your ramp should end?');
% End=1;
Delay=input('What is your Dealy?');
% Delay=0;
% decay=input('What is ypur decay rate? Should be 1 or -1');
decay=input('What is your decay? enter 1 for positive decay & -1 for negetive decay');
def=0-Min;% To plot the function on it's accurate place on x axis
Def=Min-Min;
% Stepwidth=input('How much is the step width?');
% Stepwidth=1;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
s=1/a;
d=0;% number of sequence
T0=Max-Min;
L=T0/2;
As=length(n);
ft=zeros(1,As);
% omaga=2*pi/interval;

%% find bn
omaga=2*pi/T0;

for x=1:1:length(n)
    
    if x>=(def+Delay)/a && x<=(def+End+Delay)/a% plot the ramp function
        y(x)=decay*n(x-(Delay/a));
    else
       y(x)=0;
    end
        if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)% plot the sin function for finding the value of bn
        
        z(x)=sin(d*pi*x*a/(L));
    else
        z(x)=0;
        end
        
        if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)% plot the cos function for finding the value of an
        
        w(x)=cos(d*pi*x*a/(L));
    else
        w(x)=0;
        end
end    
%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)
      ft(x)=A0;
    
else
    ft(x)=0;
end

end
%% Start the loop from n=1
for D=1:1:100000
for x=1:1:length(n) 
    if x>=(def+Delay)/a && x<=(def+End+Delay)/a
        y(x)=decay*n(x-(Delay/a));
    else
       y(x)=0;
    end
   
        if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)
        
        z(x)=sin(D*pi*x*a/(L));
    else
        z(x)=0;
        end
        
        if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)
        
        w(x)=cos(D*pi*x*a/(L));
    else
        w(x)=0;
        end
end
%% Value of A0;

we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)
      ft(x)=ft(x)+(An*cos(D*pi*x*a/(L))+Bn*sin(D*pi*x*a/(L)));
    
else
    ft(x)=0;
end

end
for x=1:1:length(n)

if x>=(Def+Delay)/a && x<=(def+Delay)/a+(Max/a)
      ftt(x)=(An*cos(D*pi*x*a/(L))+Bn*sin(D*pi*x*a/(L)));
    
else
    ft(x)=0;
end

end
RMSE=(sqrt(mean((y-ft).^2)))*100;
if RMSE<=5
    break
end

subplot(3,1,1)
  plot(n,y,'r','linewidth',2) 
  axis auto
  grid on
  xlabel('Time(in second)')
  ylabel('Amplitude')
  title('Input Ramp Function F(t)')
subplot(3,1,2)
 m=Min:(Max-Min)/Sample:Max;
plot(m,ftt,'linewidth',2)
axis auto
grid on
xlabel('Time (in second)');
ylabel('Amplitude');
  title(['Summing Signal at n=',num2str(D), ', A0=',num2str(A0),', A',num2str(D),'=',num2str(An),', B',num2str(D),'=',num2str(Bn)]);
subplot(3,1,3)
plot(n,ft,'b','linewidth',2)
axis auto
grid on
xlabel('Time(in second)')
  ylabel('Amplitude')
  title(['Output Ramp Function SF(t)',', Error in the plot is: ',num2str(RMSE)])
  pause(0.1)


end


case 3
    clc;
    close all;
e=menu('Select one Sinosudial Function','Sin Function','Cos Function');
switch e
    case 1
     clc;
        close all;
    
disp('This is Sin function:')
Min=-10;% Initialise the lower boundry
Max=10;% intialise the upper biundry of the plot
Sample=10000;% initialise the sample size
a=Sample/(Max-Min);
delay='How much Delay do you need?';
Delay=input(delay);
% Delay=-1;
def=Delay-Min;
Def=Min-Min;
Stepwidth=input('How much is the step width?');
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
Amp=input('What is your Amplitude?');
s=1/a;
% w=(0-Min)*Sample;
% Stepwidth=5;
d=0;% number of sequence
Frequency=input('What is your frequency?');
Phase=input('What is the phase? if your phase is 90 degree enter pi/2');
T0=Max-Min;
L=T0/2;;
As=length(n);
ft=zeros(1,As);
%% Value of A0;
%% find bn
omaga=2*pi/T0;
for x=1:1:length(n)
    if x>=(def)*a && x<=(def+Stepwidth)*a
          y(x)=Amp*sin(2*pi*x*Frequency/(a)+Phase);
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(d*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(d*pi*x/(a*L));
    else
        w(x)=0;
        end
end

%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=A0;  
else
    ft(x)=0;
end
end
for D=1:1:100000
for x=1:1:length(n)
    
    if x>=(def)*a && x<=(def+Stepwidth)*a
           y(x)=Amp*sin(2*pi*x*Frequency/(a)+Phase);
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(D*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(D*pi*x/(a*L));
    else
        w(x)=0;
        end
end
%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=ft(x)+(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ft(x)=0;
end

end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ftt(x)=(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ftt(x)=0;
end

end
RMSE=(sqrt(mean((y-ft).^2)))*100;
% asd=(wee-we)*100
if RMSE<=5
    break
end
subplot(3,1,1)
  plot(n,y,'r','linewidth',2) 
  axis auto
  grid on
  xlabel('Time(in second)')
  ylabel('Amplitude')
  title('Input Sin Function F(t)')
subplot(3,1,2)
 m=Min:(Max-Min)/Sample:Max;
plot(m,ftt,'linewidth',2)
axis auto
grid on
xlabel('Time (in second)');
ylabel('Amplitude');
 title(['Summing Signal at n=',num2str(D), ', A0=',num2str(A0),', A',num2str(D),'=',num2str(An),', B',num2str(D),'=',num2str(Bn)]);
subplot(3,1,3)
plot(n,ft,'b','linewidth',2)
axis auto
grid on
xlabel('Time(in second)')
  ylabel('Amplitude')
  title(['Output Sin Function SF(t)',', Error in the plot is: ',num2str(RMSE)])
  pause(0.1)
end

    case 2
       
                clc;
        close all;
    
disp('This is cos function:')
Min=-10;% Initialise the lower boundry
Max=10;% intialise the upper biundry of the plot
Sample=10000;% initialise the sample size
a=Sample/(Max-Min);
delay='How much Delay do you need?';
Delay=input(delay);
% Delay=-1;
def=Delay-Min;
Def=Min-Min;
Stepwidth=input('How much is the step width?');
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
Amp=input('What is your Amplitude?');
s=1/a;
% w=(0-Min)*Sample;
% Stepwidth=5;
d=0;% number of sequence
Frequency=input('What is your frequency?');
Phase=input('What is the phase? if your phase is 90 degree enter pi/2');
T0=Max-Min;
L=T0/2;;
As=length(n);
ft=zeros(1,As);
%% Value of A0;
%% find bn
omaga=2*pi/T0;
for x=1:1:length(n)
    if x>=(def)*a && x<=(def+Stepwidth)*a
          y(x)=Amp*cos(2*pi*x*Frequency/(a)+Phase);
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(d*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(d*pi*x/(a*L));
    else
        w(x)=0;
        end
end

%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=A0;  
else
    ft(x)=0;
end
end
for D=1:1:100000
for x=1:1:length(n)
    
    if x>=(def)*a && x<=(def+Stepwidth)*a
           y(x)=Amp*cos(2*pi*x*Frequency/(a)+Phase);
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(D*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(D*pi*x/(a*L));
    else
        w(x)=0;
        end
end
%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=ft(x)+(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ft(x)=0;
end

end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ftt(x)=(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ftt(x)=0;
end

end
RMSE=(sqrt(mean((y-ft).^2)))*100;
% asd=(wee-we)*100
if RMSE<=5
    break
end
subplot(3,1,1)
  plot(n,y,'r','linewidth',2) 
  axis auto
  grid on
  xlabel('Time(in second)')
  ylabel('Amplitude')
  title('Input Cos Function F(t)')
subplot(3,1,2)
 m=Min:(Max-Min)/Sample:Max;
plot(m,ftt,'linewidth',2)
axis auto
grid on
xlabel('Time (in second)');
ylabel('Amplitude');
 title(['Summing Signal at n=',num2str(D), ', A0=',num2str(A0),', A',num2str(D),'=',num2str(An),', B',num2str(D),'=',num2str(Bn)]);
subplot(3,1,3)
plot(n,ft,'b','linewidth',2)
axis auto
grid on
xlabel('Time(in second)')
  ylabel('Amplitude')
  title(['Output Cos Function SF(t)',', Error in the plot is: ',num2str(RMSE)])
  pause(0.1)
end
end

case 4
    clc;
    clear all;
    close all;
disp('This is Exponential function:')
% min='what is your min range?';
% Min=input(min);
Min=-10;
% max='What is your max range?';
% Max=input(max);
Max=10;
% sample='how many sample you want?';
% Sample=input(sample);
Sample=10000;
a=Sample/(Max-Min);
% delay='How much time shift do you need?';
Delay=input('What is your delay?');
% Delay=-1;
def=Delay-Min;
Def=Min-Min;
Stepwidth=input('From where your graph should stop?');
% Stepwidth=1;
n=Min:(Max-Min)/Sample:Max;% Scalling of the x axis according to sample
% Amp=1;
Amp=input('What is the Amplitude?');
s=1/a;
d=0;% number of sequence
T0=Max-Min;
L=T0/2;
As=length(n);
ft=zeros(1,As);
decay=input('What is your decay? 1 for positive decay & -1 for negetive decay');
%% Value of A0;

%% find bn
omaga=2*pi/T0;
% decay=1;
for x=1:1:length(n)
    
    if x>=(def)*a && x<=(def+Stepwidth)*a

            y(x)=Amp*exp(decay*x/a);
    else
       y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(d*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(d*pi*x/(a*L));
    else
        w(x)=0;
        end
end
we=y.*z;
yu=w.*y;
A0=trapz(n,y)/T0;
An=trapz(n,yu)/L;
Bn=trapz(n,we)/L;
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=A0; 
else
    ft(x)=0;
end
end
for D=1:1:700000
for x=1:1:length(n)
    if x>=(def)*a && x<=(def+Stepwidth)*a
            y(x)=Amp*exp(decay*x/a);
    else
           y(x)=0;
    end
        if x>=Def*a && x<=def*a+(Max*a)
        
        z(x)=sin(D*pi*x/(a*L));
    else
        z(x)=0;
        end
        
        if x>=Def*a && x<=def*a+(Max*a)
        
        w(x)=cos(D*pi*x/(a*L));
    else
        w(x)=0;
        end
end
%% Value of A0;
we=y.*z;% Multiplication of the signal f(t) & sin function to find the value of Bn
yu=w.*y;% Multiplication of the signal f(t) & Cos function to find the value of An
omaga=2*pi/L;% Find the value of Omaga
A0=trapz(n,y)/T0;% equation for A0
An=(trapz(n,yu)/L);% find the area under the curve basically find the value of An
Bn=(trapz(n,we)/L);% find the area under the curve basically find the value of Bn
if An<0.0001 && An>-0.0001
    An=0;
else
    An=An;
end
    if Bn<0.0001 && Bn>-0.0001
        Bn=0;
    else
        Bn=Bn;
    end
    An;
    Bn;
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ft(x)=ft(x)+(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ft(x)=0;
end
end
for x=1:1:length(n)

if x>=Def*a && x<=def*a+(Max*a)
      ftt(x)=(An*cos(D*pi*x/(a*L))+Bn*sin(D*pi*x/(a*L)));
    
else
    ftt(x)=0;
end
end
RMSE=(sqrt(mean((y-ft).^2)))*100;
% asd=(wee-we)*100
if RMSE<=5
    break
end
subplot(3,1,1)
  plot(n,y,'r','linewidth',2) 
  axis auto
  grid on
  xlabel('Time(in second)')
  ylabel('Amplitude')
  title('Input Exponential Function F(t)')
subplot(3,1,2)
 m=Min:(Max-Min)/Sample:Max;
plot(m,ftt,'linewidth',2)
axis auto
grid on
xlabel('Time (in second)');
ylabel('Amplitude');
 title(['Summing Signal at n=',num2str(D), ', A0=',num2str(A0),', A',num2str(D),'=',num2str(An),', B',num2str(D),'=',num2str(Bn)]);
subplot(3,1,3)
plot(n,ft,'b','linewidth',2)
axis auto
grid on
xlabel('Time(in second)')
  ylabel('Amplitude')
  title(['Output Exponential Function SF(t)',', Error in the plot is: ',num2str(RMSE)])
  pause(0.1)
end
end