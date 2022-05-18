%%
%PlayBack calls
sigamask=bandpass(sig1mask,[35000 70000],fs);
t1mask = [1:1:numel(sigamask)]*1/fs;

siga1call=sigamask(358000:359100); %the strongest call in Sophi file- willd o cross corerrelation on it.
[cross_corr, lags] = xcorr(sigamask,siga1call);
tao = lags*1/fs;

[pks,locs]=findpeaks(cross_corr,'MinPeakDistance',1.2*fs,'MinPeakHeight',max(cross_corr)/24);%1.2*fs);
tao2 = locs*(1/fs)-t1mask(end);
taodiff=diff(tao2);

figure
plot(t1mask,sig1mask)
title('playback calls')
hold on
findpeaks(cross_corr,'MinPeakDistance',1.2*fs,'MinPeakHeight',max(cross_corr)/24);

%%
%No-Noise Trial
sigb=bandpass(sig2shNN,[35000 70000],fs);
t2shani = [1:1:numel(sigb)]*1/fs;

[cross_corrShani, lagsShani] = xcorr(sigb,siga1call);
taoShani = lagsShani*1/fs;

[pkshani,locshani]=findpeaks(cross_corrShani,'MinPeakDistance',0.002*fs, 'MinPeakHeight',max(cross_corr)/48);
tao2shani = locshani*(1/fs)-t2shani(end);
diffShani=diff(tao2shani);
meanshani=mean(diffShani);
ipiNoNoise=meanshani;
%%
%Trial with conspecific's echolocation calls
sigbCon=bandpass(sig2sh,[35000 70000],fs);
t2shani = [1:1:numel(sigbCon)]*1/fs;

[cross_corrShani, lagsShani] = xcorr(sigbCon,siga1call);
taoShani = lagsShani*1/fs;

[pkshani,locshani]=findpeaks(cross_corrShani,'MinPeakDistance',0.002*fs, 'MinPeakHeight',max(cross_corr)/48);
tao2shani = locshani*(1/fs)-t2shani(end);
diffShani=diff(tao2shani);
meanshani=mean(diffShani);

taodiff=diff(tao2);
csumdiff=cumsum(taodiff);
dtime = [];
l=[];
timesinShani=[];
taobeforecall=cumsum(taodiff(end-2:-1:1));
taobeforecall=taobeforecall(end:-1:1);
for ii=1:length(tao2shani)-length(taodiff)-1
    if tao2shani(ii)+2<=tao2shani(end)
        smallertao=tao2shani(tao2shani<=tao2shani(ii)+2);

        %I will use ismembertol for understanding
        %where the playback calls are inside the signal.
        %I will compare the time interval of the
        %playback calls and the peaks
        %(from the current call, and 2 seconds ahead-
        %because the playback signal is only 1.2 sec) in the signal.
        %I will throw every call depicted in the signal.
%##########
        aa=ismembertol([tao2shani(ii:find(tao2shani==smallertao(end)))],[tao2shani(ii);tao2shani(ii)+csumdiff],0.003,'DataScale', 1);
        l=[l,sum(aa)];
        if sum(aa)<5
            continue
        else

            dtime= [dtime; tao2shani(ii)];
            timesinShani=[timesinShani; tao2shani(ii)-taobeforecall', tao2shani(ii),tao2shani(ii)+cumsum(taodiff(13:14))']
        end
    else
        aa=ismembertol([tao2shani(ii:find(tao2shani==smallertao(end)))],[tao2shani(ii);tao2shani(ii)+csumdiff],0.003,'DataScale', 1);

        if sum(aa)<5
            continue
        else
            dtime= [dtime; tao2shani(ii)];
            timesinShani=[timesinShani; tao2shani(ii)-taobeforecall', tao2shani(ii),tao2shani(ii)+cumsum(taodiff(13:14))'];
        end
    end
end
sigblesss=sigbCon;
[lia,locsbless]=ismembertol(timesinShani,t2shani) ;
locsbless = reshape(locsbless,[],1);

sortlocs=sort(locsbless);
sortlocs=sortlocs(sortlocs>0);
sigblesss(sortlocs)=0;
lo=sigbCon(sortlocs)~=sigblesss(sortlocs);
sumlo=sum(lo);
ipiBistatic=[meanshani];

figure()
plot(t2shani,sigbCon)
hold on
plot(tao2shani,pkshani/(max(pks)*2),'o')
title('peaks of cross-correlation on the signal')
xlabel('time(sec)')
ylabel('power')
legend('signal','peaks')

%%
%FFTT GRAPH
sig2call=sigbCon(3555800:3556900);
siga1call=sigamask(358000:359100);

f1=fs*(0:(fs/4)-1)/(fs/2); % frequency axis
Y=fft(sig2call,fs/2); 
Pyy=Y.*conj(Y)/(fs/2); % Power spectrum
figure % Plot result
plot(f1,Pyy(1:(fs/4)));
f1(1)=0;
title('Powerspectrum')
xlabel('Frequency (Hz)')
ylabel('Power (mV2)')
hold on

f=fs*(0:(fs/4)-1)/(fs/2); % frequency axis
Y2=fft(siga1call,fs/2);
Pyy=Y2.*conj(Y2)/(fs/2); % Power spectrum
plot(f,Pyy(1:(fs/4)));
f(1)=0;
legend('Shani calls', 'Playback calls')

%%
%statistics

[hf,pf] = vartest2(ipiNoNoise,ipiBistatic) 
[h,p] = ttest2(ipiNoNoise,ipiBistatic,'Vartype','unequal') 

figure()
bar(ipiNoNoise,'r')
hold on
bar(ipiBistatic,'b')
title('Bar Plot Of Mean IPI')
xlabel('Trial')
ylabel('Mean IPI')
legend('IPI-NoNoise', 'IPI-Bistatic')

%%
%all the statistics:
figure()
bar(ipi1NoNoise,'r')
hold on
bar(ipiBiatatic,'b')
title('Bar Plot Of Mean IPI')
xlabel('Trial')
ylabel('Mean IPI')
legend('IPI-NoNoise', 'IPI-Bistatic')