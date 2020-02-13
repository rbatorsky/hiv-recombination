% Updated Jan 2011
% plots 1 fig and 1 txt file 
% changed recomb model (1/11)

% for Fig. 2 of observables vs. r and s0
% recomb_0('constant',0.05,0.1,1,500,1000,10,1500,1)

function recomb_0(distribution_s,r,s0,a,L,N,M,tf,run)
homedir = '/cluster/shared/rbator01/recomb_model_test/new_model/';
%homedir = '~/Desktop/';
filename = sprintf('%s_r=%g_s0=%g_%g_%g_%g_%g_%g',...
    distribution_s,r,s0,a, M,L,N,run);
%diary on
%diary(strcat(homedir,filename,'_diary.txt'))

% random seed, so can average over runs.

 RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

f0=0;                         % preset initial frequency of alleles
muesc =0;                     % standing variation f0(s): deleterious mutation rate during escape
mu= 3e-5;                     % recurrent beneficial mutation rate 
tesc=200;                     % time interval of escape
%M=10;                        % crossovers per genome
fmin=0.04;                    % cutoff  for average sample of 25 seq. 
T=0:tf';                      % times
Lsegment=round(L/9);          % Sampling segment, as in data

% Distributions of s

switch distribution_s

    case 'power' 
        % 1: a*s0^a/(s+s0)^(a+1) 
        s=s0*(rand(1,L).^(-1/a)-1); 
        ts=sprintf('1/s^{a+1},s>s0. N=%g,r=%g,L=%g,s0=%g,\n a=%g,f0=%g,fmin=%g,T=%g,M=%g',N,r,L,s0,a,f0,fmin,max(T),M);

    case 'lognormal'
        % 1/s*exp[-log^2(s/s0)/2/a^2]
        s=s0*exp(a*randn(1,L));
        ts=sprintf('1/s*exp[-log^2(s/s0)/2/a^2]\n N=%g,r=%g,L=%g,s0=%g,a=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n fmin=%g,T=%g,M=%g',N,r,L,s0,a,f0, muesc,mu,tesc,fmin,max(T),M);
    
    case 'exponential'
        % exp(-s/s0) %% actually density is, g(s) = (1/s0)*exp(-s/s0)
        s=-s0*log(rand(1,L));
        ts=sprintf('exp(-s/s0)\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n fmin=%g,T=%g,M=%g',N,r,L,s0,f0, muesc,mu,tesc,fmin,max(T),M);%
    
    case 'constant'
        s = s0*ones(1,L);
        ts=sprintf('constant\n N=%g,r=%g,L=%g,s0=%g\n f0=%g, muesc=%g,mu=%g,tesc=%g,\n fmin=%g,T=%g,M=%g',N,r,L,s0,f0, muesc,mu,tesc,fmin,max(T),M);%
    
end


ns=10;                                                  % number of log s bins for plotting
smax=max(s); smin=min(s);
sbound=exp(log(smin)+(0:ns)*(log(smax)-log(smin))/ns);  % bin boundaries
disp(ts)
disp(sprintf('exp(sum(s))=%g, smax=%g',exp(sum(s)),smax))  % max possible fitness increase

% Initial settings
figure(1);subplot(2,2,1);clf;subplot(2,2,2);clf;subplot(2,2,3);clf;subplot(2,2,4);clf;
figure(2);subplot(2,2,1);clf;subplot(2,2,2);clf;subplot(2,2,3);clf;subplot(2,2,4);clf;
fsite=zeros(length(T),L); 
ndivbin = zeros(length(T),ns); favbin=zeros(length(T),ns); davbin=zeros(length(T),ns); 
Knew=zeros(N,L);

% x should be a factor of tf in order to give x
% times plotted in figure for wave and dav. Take tf = 1500 and x = 5.

tint=tf/5; 


% for L/9 (first 9th of genome to mimic pol data)
zs=zeros(length(T),1);
fav=zs; dav=zs;  obsdiv=zs; vdiv = zs;
favseg=zs; davseg=zs; obsdivseg=zs; vdivseg = zs; 
LD4 = [];f3haplow = [];fitnessgain = [];
LD4_tf = 0; f3haplow_tf = 0;
% Matrix K: Each row is a genome, sequence of 0 and 1
%% Initial population: 
% Case 1: randomly distributed good alleles with fixed frequency
if f0~=0
    K=(rand(N,L) < f0); 
% Case 2: dynamic distribution based on the antigenic escape compensation model: 
% Alleles acumulate during tesc as deleterious, then change the sign of s. 
% Model valid if 1-site deterministic: f0(s) > 1/Ns <=> mu*N > 1 for s > 1/tesc)
elseif muesc~=0
	f0=muesc*s.^(-1).*(1-exp(-s*tesc));       % row of initial frequencies 
	K=(rand(N,L) < ones(N,1)*f0);
% Case 3: no alleles at start
else
    K=zeros(N,L);
end

% Evolution starts...
for t=T
    % beneficial mutation
    K = K | sprand(N,L,mu);
    % recombination of some pairs with parent replacement
    npairs=round(r*N/2);
    ii=ceil(rand(npairs,2)*N); i1=ii(:,1); i2=ii(:,2);  % 2 columns of random indices
    for i=1:npairs
        % generating random 1 or 0 for each site, with prob. M/L and 1-M/L,
        % respectively; even/odd xx shows site segments copied from 1st parent,
        % 2nd, 1st, 2nd etc
        xx=cumsum(rand(1,L) < M/L); 
        first=(round(xx/2)==xx/2);    % sites copied from 1st parent
        prog1=K(i1(i),:).*first+K(i2(i),:).*(1-first);    % progeny sequence
        %prog2=K(i2(i),:).*first+K(i1(i),:).*(1-first);    % throw away complementary progeny
        K(i1(i),:)=prog1; %K(i2(i),:)=prog2;               % only one parent is replaced (1/11)
    end
    
    % sampling and selection
    % column of N log-fitnesses 
    % state with 0s only has w=0 (fitness 1) by definition
    w=K*(s'); 
    % Random drift: divide in 2 or die 
    % Same as Poisson, Var[n]=E[n], Neff=N for drift/tree; not for few copies of an allele
    nprogav=min(2,exp(w)/sum(exp(w)));      % average progeny numbers <= 2
    nprogav=nprogav/mean(nprogav);          % renormalizing
    nprog=rand(N,1) < nprogav/2;            % survive or die?
    % balancing survivals and deaths to keep population constant
    disbal=1;
    while disbal ~=0
        isur=find(nprog); idied=find(1-nprog);
        disbal=length(isur) - N/2;
        if disbal > 0
            iflip=isur(ceil(rand(1,disbal)*length(isur)));
            nprog(iflip)=zeros(1,disbal); % replacing extra 1 by 0
        elseif disbal < 0
            iflip=idied(ceil(rand(1,-disbal)*length(idied)));
            nprog(iflip)=ones(1,-disbal);
        end
    end
    nprog=nprog*2;
    % updating population
    is=[0;cumsum(nprog(1:(N-1)))];
    for i=1:N
        if nprog(i)
            Knew(is(i)+1:is(i)+nprog(i),:)=ones(nprog(i),1)*K(i,:);
        end
    end
    K=Knew; 
    sK=size(Knew);sK=sK(1);
    if sK  ~= N, disp('N not conserved'),return;end
    
    % Evolution step ended
    
    % Memorizing various values
    % Memorizing average allele frequency for each site
    meanK=mean(K); 
    fsite(t+1,:)=meanK;
    % ... and for each s bin
    for k=1:ns
        favbin(t+1,k)=mean(meanK(s > sbound(k) & s <= sbound(k+1)));
        ii=find(s > sbound(k) & s <= sbound(k+1) & meanK > fmin & meanK < 1-fmin);
        ndivbin(t+1,k) = length(ii);
        davbin(t+1,k)=mean(2*meanK(ii).*(1-meanK(ii)));     % average diversity for bin over observ. sites
    end
    
    
    % record observables, for whole genome 
    ii=find(fsite(t+1,:) > fmin & fsite(t+1,:) < 1-fmin);   % observably diverse sites
    iii = find(fsite(t+1,:) > 0.25 & fsite(t+1,:) < 0.75);
    obsdiv(t+1)=length(ii)/L;                               % their frequency for each time
    vdiv(t+1) = length(iii)/L;
    fav(t+1)=mean(fsite(t+1,ii));                           % their average-over-site frequency
    dav(t+1)=mean(2*fsite(t+1,ii).*(1-fsite(t+1,ii)));      % same for diversity
    
    % record observables, for L/9
    ii=find(fsite(t+1,1:Lsegment) > fmin & fsite(t+1,1:Lsegment) < 1-fmin);   % observably diverse sites
    iii = find(fsite(t+1,1:Lsegment) > 0.25 & fsite(t+1,1:Lsegment) < 0.75);
    obsdivseg(t+1)=length(ii)/Lsegment;                               % their frequency for each time
    vdivseg(t+1) = length(iii)/Lsegment;    
    favseg(t+1)=mean(fsite(t+1,ii));                           % their average-over-site frequency
    davseg(t+1)=mean(2*fsite(t+1,ii).*(1-fsite(t+1,ii)));      % same for diversity
    

    %%
    % Plotting results at 11 time points
    %%
    
    col='rgbymckrgbymckrgbymckrgbymck';
    if tint*round(t/tint)==t
        c=col(round(t/tint)+1);
        
        figure(1)
        subplot(2,2,1)
        hold on
        [nn,xx]=hist(exp(w));                % histogram of fitness among genomes
        plot(xx,nn,c)
        %{
        subplot(2,2,2)
        hold on
        [nn,xx]=hist(sum(K,2));         % histogram of total k among genomes
        plot(xx,nn,c)
        %}
        subplot(2,2,4)
        % Average diversity for each bin for observably diverse sites
        bincenter=(sbound(1:end-1)+sbound(2:end))/2;
        semilogx(bincenter,davbin(t+1,:),c)        
        hold on
    
        figure(2)
        subplot(1,2,1)
        bincenter=(sbound(1:end-1)+sbound(2:end))/2;
        semilogx(bincenter,favbin(t+1,:),c)   
        hold on
        
        subplot(1,2,2)
        bincenter=(sbound(1:end-1)+sbound(2:end))/2;
        semilogx(bincenter,ndivbin(t+1,:),c)
        hold on
        
    end
        % Clone and haplotype structure at same times except t=0
    if t~=0 && tint*round(t/tint)==t 

            % Size-ranked haplotype frequencies for v.div.site pairs at this showtime, same as above for tstarthap
            idiv=find(meanK(1:Lsegment) > 0.25 & meanK(1:Lsegment) < 0.75);  % Selecting very diverse sites
            ss=length(idiv);
        if ss >= 3   
                nphap=ss*(ss-1)/2;         % total v.d.pairs at this time
                
                 % Matrix of locations of v.div. pairs
                 % changed to use nchoosek 8/7
            
            ips = nchoosek(idiv,2);
            fhap=zeros(nphap,4); 
            fhapnorm=zeros(nphap,4);
            n3hap = 0;n3haplow = 0;
            
            %d = zeros(nphap,2); %[distance, LD] for binning plot

            % loop over pairs
            for i=1:nphap
                 
                 f00=mean(all(K(:,ips(i,:))==ones(N,1)*[0 0],2));
                 f01=mean(all(K(:,ips(i,:))==ones(N,1)*[0 1],2));
                 f10=mean(all(K(:,ips(i,:))==ones(N,1)*[1 0],2));
                 f11=mean(all(K(:,ips(i,:))==ones(N,1)*[1 1],2));
                 [fhap(i,:),ii]=sort([f00 f01 f10 f11],'descend');    % ranked hap. frequency
                 f1=mean(K(:,ips(i,1)));         % 1-site freq. for site 1 of the pair
                 f2=mean(K(:,ips(i,2)));         % same for site 2
            % order of denominators now correct, fixed on July 18th
                 denom=[(1-f1)*(1-f2) (1-f1)*f2 f1*(1-f2) f1*f2];    % expected w/o LD
                 denom=denom(ii);
                 fhapnorm(i,:)=fhap(i,:)./denom; % norm. ranked hap. freq.
                 
                 
                 hapsmissing =   (f00 == 0) + (f01 == 0) +  (f10 == 0) +  (f11 == 0);
                 hapslow =   (f00 < fmin) + (f01 < fmin) +  (f10 < fmin) +  (f11 <fmin);               

                 n3hap = (hapsmissing  ==  1) + n3hap; 
                 n3haplow = (hapslow  ==  1) + n3haplow;
            end    
                figure(1)
                subplot(2,2,2);hold on
 
                plot(1:4,mean(fhap),[c 'o'],1:4,mean(fhap),c,...
                1:4,mean(fhapnorm),[c 'x'],1:4,mean(fhapnorm),c,...
                [1 2 3 4;1 2 3 4],[mean(fhapnorm)-std(fhapnorm);mean(fhapnorm)+std(fhapnorm)],[c 'x'])
                                  
            
            
                f3haplow = [f3haplow n3haplow/nphap];
                LD4 = [LD4 mean(fhapnorm(:,4))];

            if t == tf
                LD4_tf = mean(fhapnorm(:,4));
                f3haplow_tf = n3haplow/nphap;
            end
            
        else  
                disp(sprintf('no v. div. pairs at t=%g',t))
        end
    end
        
    if t == tf
        fitnessgain = mean(exp(w));
    end
    

        %end %  tf    
    
end %Evolution ended



%Final plots
figure(1)

subplot(2,2,1)
hold off
%axi=axis; axi(1:2)=[1 exp(sum(s))]; axis(axi);
xlabel('fitness exp(w) ');
title(ts)   % title  with parameter values

subplot(2,2,3)
plot(T,fav,'b',T,dav,'g',T,obsdiv,'m',T, vdiv, 'k')
title({'whole genome';'<f> (b) and <2f(1-f)> (g) for obs.div. sites;';' site fract.: obs.div. (m), v.div. (k), not fixed (r)'})
xlabel('t');ylabel('fav, dav')
axi=axis;axi(3:4)=[0 0.5];axis(axi);


if smin ~= smax
    
            subplot(2,2,4)
            hold off
            xlabel('s'); ylabel('dav_s');
            %axi=axis;axi(1:2)=[smin smax];axis(axi);
            title(sprintf('Avrg. (over observably diverse\n sites) diversity of %g bins',ns))
            

            
end

figure(1)
saveas(gcf,strcat(homedir,filename,'.fig'));
figure(2)
saveas(gcf,strcat(homedir,filename,'_2.fig'));

obsdivseg = obsdivseg';
vdivseg = vdivseg';
davseg = davseg';



%recording info
switch distribution_s
    
    case 'exponential'
        distribution = 1 ;
    case 'lognormal'
        distribution = 2;
    case 'constant'
        distribution = 3;
    case 'power'
        distribution = 4;
end

trecord = 1500;
info1500 = sprintf('%g \n',distribution, r, s0, N, LD4_tf,f3haplow_tf,obsdivseg(trecord),vdivseg(trecord),davseg(trecord),fitnessgain,M);
dlmwrite(strcat(homedir,'txt/',filename,'.txt'),info1500,'');

end




