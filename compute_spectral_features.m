%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Jesús Monge-Álvarez/ Patricia Amado Caballero
% Email: patricia.amado@uva.es
% Date: 2025-08-27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RP_feat,SpecBand_feat,SpecCent_feat,SpecCrestFac_feat,...
    SpecEn_feat,SpecFlat_feat,SpecFlux_feat,SpecKurt_feat,...
    SpecRenyiEn_feat,SpecRolloff_feat,SpecSkew_feat]=...
    calcula_features_espectrales(spectrogram,fs,bandLim,t,time_lim)

% Parameters:
%  - spectrogram: one-sided PSD, frequencies in rows, time in columns
%  - bandLim: frequency band limits to analyze (numBands x 2 matrix)
%  - time_lim: time segment limits (numT x 2 matrix)
%
% Returns:
%  - Set of spectral features of size numT x numBands
%    (except SpecEn, which is of size numT x 1)


frames=spectrogram';
nfft=size(frames,2);
f = linspace(0,fs/2,nfft);

numBands = size(bandLim,1);
indBands = NaN(numBands,numel(f));

for IDband = 1:1:numBands
    cfmin = bandLim(IDband,1);
    cfmax = bandLim(IDband,2);
    
    tmp = (f>=cfmin) & (f<cfmax);
    
    indBands(IDband,:) = tmp;
    
    clear('cfmin','cfmax','tmp');
end % end for each freq. band


numT=size(time_lim,1);
indT = NaN(numT,numel(t));


for IDT = 1:1:numT
    cTmin = time_lim(IDT,1);
    cTmax = time_lim(IDT,2);
    
    tmp = (t>=cTmin) & (t<cTmax);
    
    indT(IDT,:) = tmp;
    
    clear('cTmin','cTmax','tmp');
end % end for each freq. band


%%
numFra = size(frames,1);
Ratiof50vsf90 = NaN(numFra,numBands);

for IDfra = 1:1:numFra
    cframe = frames(IDfra,:);
    if IDfra>1
        cframe2=frames(IDfra-1,:);
    end

    
    % Welch's PSD estimation:
   % wPSD = pwelch(cframe,WelchWlen,WelchOverlap,nfft,fs);
   % wPSD = wPSD(:).';
    aux = NaN(1,numBands);

    for IDband = 1:1:numBands
        cindx = logical(indBands(IDband,:));
        tmp = cframe(cindx); 
        cf = f(cindx);

        % Ratiof50vsf90 computation:
        totalEne = sum(tmp);
        counter = 1;
        sumEne = 0;
        
        limit_50 = 0.50*totalEne;
        limit_90 = 0.90*totalEne;
        
        if limit_50>0
            a=1;
        end
        
        while sumEne<limit_50
            sumEne = sumEne + tmp(counter);
            
            counter = counter + 1;
        end
        if counter==1
            counter=2;
        end
        f50 = cf(counter - 1);
        
        while sumEne<limit_90
            sumEne = sumEne + tmp(counter);
            
            counter = counter + 1;
        end
          if counter==1
            counter=2;
        end
        f90 = cf(counter - 1);
        
        Ratiof50vsf90(IDfra,IDband) = f90/f50;
        clear('totalEne','sumEne','counter', ...
            'limit_50','limit_90','f50','f90');
        %SpecBand calculation
        Num_SpecCent = sum(cf.*tmp);
        Den_SpecCent = sum(tmp);
        SpecCent = Num_SpecCent/Den_SpecCent;
        
        Term1 = (cf - SpecCent).^2;
        Num = sum(Term1 .* tmp);
        Den = Den_SpecCent;
        
        SpecBand(IDfra,IDband) = Num/Den;
        clear('Num_SpecCent','Den_SpecCent','SpecCent', ...
        'Term1','Num','Den');

         %SpecCent calculation
         Num = sum(cf.*tmp);
         Den = sum(tmp);
         SpecCent(IDfra,IDband) = Num/Den;
         clear('Num','Den');
        
         
         %SpecCrestFact calculation
         cte = max(cf) - min(cf) + 1;
         Num = max(tmp);
         Den = (1/cte) * sum(tmp);
         SpecCrestFac(IDfra,IDband) = Num/Den;
         clear('Num','Den');

         %RP calculation
         TotalEne=sum(cframe);       
         cBandEne = sum(tmp);
         RP(IDfra,IDband) = cBandEne/TotalEne;

         clear('cBandEne')

     % SpecFlat computation:
          
        Term1 = mean(log(tmp));
        Num = exp(Term1);
        Den = mean(tmp);
        SpecFlat(IDfra,IDband) = Num/Den;
        clear('Term1','Num','Den');
  
    %SpecKurt computation
        Term1 = 10*log10(tmp);
        Term2 = mean(Term1);
        
        Num = mean((Term1 - Term2).^4);
        Den = (std(Term1,1).^4);
        
        SpecKurt(IDfra,IDband) = Num/Den;
         clear('Num','Den');

      %SpecSkew computation
              
        Num = mean((Term1 - Term2).^3);
        Den = (std(Term1,1).^3);
        
        SpecSkew(IDfra,IDband) = Num/Den;
        
       clear('Term1','Term2','Num','Den');



         %SpecPeakEnd computation (we do not use it. Too many NaNs if peaks
         %are not found
        % Normalization:
        %tmp2 = tmp/sum(tmp);
        
        % Find the local maxima:
        %pks = findpeaks(tmp2);
        
      %  SpecPeakEn(IDfra,IDband) = (-1)*sum(pks .* log10(pks));
         

       %    clear('tmp2','pks');

       %SpecRenyiEn computation
        q=4; 
        cte = 1/(1-q);
        Term1 = sum(tmp.^q);
        SpecRenyiEn(IDfra,IDband) = cte * log(Term1);
        clear('q','cte','Term1')
         % SpecEn computation:
         aux(IDband) = sum(tmp)/TotalEne;
       

         %SpecRollOff computation
        TotalEne = sum(tmp);
        limit = 0.85*TotalEne;
        counter = 1;
        sumEne = 0;
        
        while sumEne<limit
            sumEne = sumEne + tmp(counter);
            
            counter = counter + 1;
        end
        if counter == 1
            counter=2;
        end
        SpecRolloff(IDfra,IDband) = cf(counter - 1);
        clear('sumEne','limit','counter','TotalEne','cf');

        %SpecFlux computation
        if IDfra==1
            SpecFlux(IDfra,IDband) = sum((tmp).^2);
        else
            tmp2 = cframe2(cindx); 
            SpecFlux(IDfra,IDband) = sum((tmp - tmp2).^2);
        end
              clear('tmp2');


          clear('cindx','tmp');
    end % end for each freq. band
    
        Term1 = aux .* log2(aux);
        SpecEn(IDfra) = (-1) * sum(Term1);
        clear('aux','Term1')

 
end % end for each frame
SpecCent_feat=NaN(numT,numBands);
SpecBand_feat=NaN(numT,numBands);
Ratiof50vsf90_feat=NaN(numT,numBands);
SpecCrestFac_feat=NaN(numT,numBands);
RP_feat=NaN(numT,numBands);
SpecFlat_feat=NaN(numT,numBands);
SpecKurt_feat=NaN(numT,numBands);
SpecSkew_feat=NaN(numT,numBands);
%SpecPeakEn_feat=NaN(numT,numBands);
SpecRenyiEn_feat=NaN(numT,numBands);
SpecRolloff_feat=NaN(numT,numBands);
SpecFlux_feat=NaN(numT,numBands);
SpecEn_feat=NaN(numT,1);



  for IDT = 1:1:numT
        cindx = logical(indT(IDT,:));
        tmp = SpecCent(cindx,:);
        SpecCent_feat(IDT,:)=mean(tmp,1);
        tmp = SpecBand(cindx,:);
        SpecBand_feat(IDT,:)=mean(tmp,1);
        tmp = Ratiof50vsf90(cindx,:);
        Ratiof50vsf90_feat(IDT,:)=mean(tmp,1);
        tmp = SpecCrestFac(cindx,:);
        SpecCrestFac_feat(IDT,:)=mean(tmp,1);
        tmp = RP(cindx,:);
        RP_feat(IDT,:)=mean(tmp,1);
         tmp = SpecFlat(cindx,:);
        SpecFlat_feat(IDT,:)=mean(tmp,1);
        tmp = SpecSkew(cindx,:);
        SpecSkew_feat(IDT,:)=mean(tmp,1);
        tmp = SpecKurt(cindx,:);
        SpecKurt_feat(IDT,:)=mean(tmp,1);
         %   tmp = SpecPeakEn(cindx,:);
     %   SpecPeakEn(IDT,:)=mean(tmp,1);
           tmp = SpecRenyiEn(cindx,:);
        SpecRenyiEn_feat(IDT,:)=mean(tmp,1);  
           tmp = SpecRolloff(cindx,:);
        SpecRolloff_feat(IDT,:)=mean(tmp,1); 
            tmp = SpecFlux(cindx,:);
        SpecFlux_feat(IDT,:)=mean(tmp,1);
     tmp = SpecEn(cindx);
        SpecEn_feat(IDT,:)=mean(tmp);
  end





    

