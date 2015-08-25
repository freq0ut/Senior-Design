close all;
clear all;
clc;

fS = 100;
tS = 1/fS;

fSce = 1;

N0 = 1024;

td1 = (2*rand()-1);
td2 = (2*rand()-1);
td3 = (2*rand()-1);
td4 = (2*rand()-1);

t = 0:tS:(N0-1)*tS;
y1 =  7+cos(2*pi*fSce*(t+td1));
y2 =  5+cos(2*pi*fSce*(t+td2));
y3 =  3+cos(2*pi*fSce*(t+td3));
y4 =  1+cos(2*pi*fSce*(t+td4));

TOA1 = 4+2*rand()-1;
TOA2 = 4+2*rand()-1;
TOA3 = 4+2*rand()-1;
TOA4 = 4+2*rand()-1;

for i=1:N0;
   if ( i < round(TOA1/tS) )
       y1(i) = 7;
   end
   
   if ( i < round(TOA2/tS) )
       y2(i) = 5;
   end
   
   if ( i < round(TOA3/tS) )
       y3(i) = 3;
   end
   
   if ( i < round(TOA4/tS) )
       y4(i) = 1;
   end
end

figure(1);
    plot(t,y1,'b')
    hold on;
    plot(t,y2,'r')
    plot(t,y3,'g')
    plot(t,y4,'m')
    ylim([-1,9]);
    legend('Chan1','Chan2','Chan3','Chan4');
    hold off;
