function err = fitERF2(start,dummy, THESIZE)
%global dummy THESIZE

%i_middle = start(1); 
bbb = start(1); 
offRr = start(2);
Cmax = start(3);

xxtemp = [1:length(dummy)];
%ppp_erf1 = erf(( xxtemp   - i_middle)/bbb)*max(dummy);
%ppp_erf1 = erf(( xxtemp   - i_middle)/bbb);
%ppp_erf1 = erf((xxtemp - round(THESIZE/2))/bbb);
ppp_erf = offRr + 1./(Cmax+exp(bbb.*xxtemp));

% whos ppp_erf dummy
% cvbvcb
err = sum(abs(sqrt(ppp_erf-dummy').^2 ));

% if offRr <0
%     err = err*10;
% end

if 0
    figure(65);clf
    plot(xxtemp,dummy,'k*')
    hold on
    plot(xxtemp,ppp_erf,'b-')
%     plot(xxtemp,ppp_erf1,'r-')
%     plot(xxtemp,ppp_erf2,'r--')
    drawnow

end



%pause(0.12)