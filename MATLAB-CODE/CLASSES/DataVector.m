classdef DataVector
    properties
        d
        Locations
        RayParams
        BackAzimuths
        InputDataParams
        DataParams
        nptsPerSeis
    end
    methods
        function obj=DataVector(InputDataParams)
            
            %Subfunctions
            function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices(InputDataParams)

                skipSta=InputDataParams.skipSta;
                KTimes= InputDataParams.KTimes;
                incangs=InputDataParams.Angles;
                TakeDifferences=InputDataParams.TakeDifferences;
                DeconvolveParentWaveform=InputDataParams.DeconvolveParentWaveform;
                obj.DataParams=InputDataParams;

                Sta1=1;
                Sta2=150;
                Stas=Sta1:skipSta:Sta2;

                if (skipSta < 0);
                    Stas=Sta2:skipSta:Sta1;
                end

                obj.nptsPerSeis=length(KTimes);
                
                nSeis=length(Stas)*length(incangs);

                %initialize matrices
                Data=zeros(nSeis,obj.nptsPerSeis);
                Locations=zeros(nSeis,1);
                RayParams=zeros(nSeis,1);
                BackAzimuths=zeros(nSeis,1);

                iseis=0;

                fprintf('Reading files... ')
                counter=0;
                for incang = incangs;

                    ModelDirectory=sprintf(...
                        'TEST_MODELS/L2_DEPTH/OUTPUT_FILES_%d-%d-%d-%d/',...
                        incang,...
                        InputDataParams.lab_amplitude,...
                        InputDataParams.lab_wavelength,...
                        InputDataParams.lab_depth);
                    
                    RfnceDirectory=sprintf(...
                        'TEST_MODELS/L2_DEPTH/OUTPUT_FILES_%d-%d-%d-%d/',...
                        incang,...
                        0.0,...
                        InputDataParams.lab_wavelength,...
                        InputDataParams.lab_depth);
                    
                    command=['cat ' ModelDirectory '/DATA/STATIONS | awk {''print $3''} > tmp.txt'];
                    system(command);
                    AllLocations=load('tmp.txt');
                    AllLocations=AllLocations/double(1000.0); % convert to km


                    for iSta = Stas;
                        iseis=iseis+1;
                        Daughter=load(sprintf([ModelDirectory 'OUTPUT_FILES/AA.S0%03d.BXP.semd'],iSta));
                        Parent  =load(sprintf([ModelDirectory 'OUTPUT_FILES/AA.S0%03d.BXS.semd'],iSta));

                        %PRef    =load(sprintf([RfnceDirectory 'OUTPUT_FILES/AA.S0%03d.BXS.semd'],iSta));

                        [AmpMax,imax]=max(abs(Parent(:,2)));
                        time=Parent(:,1) - Parent(imax,1);

                        %Parent(:,2)=Parent(:,2)/AmpMax;
                        Daughter(:,2)=Daughter(:,2)/AmpMax;

                        %plot(Parent(:,1),Parent(:,2),Daughter(:,1),Daughter(:,2))
                        %pause

                        %Subtract reference waveform
                        if (TakeDifferences)
                            DRef    =load(sprintf([RfnceDirectory 'OUTPUT_FILES/AA.S0%03d.BXP.semd'],iSta));
                            DRef(:,2) = DRef(:,2)/AmpMax;
                            Daughter(:,2)=Daughter(:,2)-DRef(:,2);
                            %Parent(:,2)  =Parent(:,2)  -PRef(:,2);
                        end

                        if (DeconvolveParentWaveform)
                           Mask=(Parent(:,1)>Parent(imax-40,1)).*(Parent(:,1)<Parent(imax+40,1));
                           P=Parent(:,2).*Mask;
                           D=Daughter(:,2);
                           TB=4;
                           NT=3;
                           dt=Parent(2,1)-Parent(1,1);
                           win_len=50;
                           Poverlap=0.99;        

                           [Time, RF_Time] = ETMTM(P(1:imax+200)',D(1:imax+200)',TB,NT,'data',dt,win_len,Poverlap);

                           %clf;
                           %subplot(2,1,1)
                           %plot(P); hold on
                           %plot(D);
                           %subplot(2,1,2)
                           %plot(Time,RF_Time);
                           %pause

                           %reset and resize daughter
                           Daughter=zeros(length(RF_Time),2);
                           Daughter(:,2)=RF_Time;
                           time=Time;

                        end

                        dtmp=interp1(time,Daughter(:,2),KTimes)';

                        Data(iseis,:)=dtmp;
                        BackAzimuths(iseis,1)=1.0;
                        RayParams(iseis,1)=incang;
                        Locations(iseis,1)=AllLocations(iSta);
                        
                        counter=counter+1;
                        fprintf('%5.2f pct\n',counter/nSeis*100.0)

                    end
                end

                fprintf('...done!\n ')
                
            end

            function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_reverse(InputDataParams)

                InputDataParams.skipSta=-InputDataParams.skipSta;

                [Data,Locations,RayParams,BackAzimuths]=build_input_matrices(InputDataParams);

                %Reflect locations about symmetry axis
                x0=1725.0;
                tmp= 2*x0 - Locations;

                Locations = tmp;

                BackAzimuths=BackAzimuths*-1.0;

            end

            function [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_both_directions(InputDataParams)
                [Data1,Locations1,RayParams1,BackAzimuths1] = build_input_matrices(InputDataParams);
                [Data2,Locations2,RayParams2,BackAzimuths2] = build_input_matrices_reverse(InputDataParams);
                Data = [Data1; Data2];
                Locations = [Locations1; Locations2];
                RayParams = [RayParams1; RayParams2];
                BackAzimuths = [BackAzimuths1; BackAzimuths2];
            
            end               
                
            function [Data,Locations,RayParams,BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,factor)

                nData=length(Locations);

                mask=rand(nData,1)>factor;

                Data=Data(mask,:);
                Locations=Locations(mask);
                RayParams=RayParams(mask);
                BackAzimuths=BackAzimuths(mask);

            end

            fprintf('Building input matrices... ')
            if InputDataParams.direction==1;
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices(InputDataParams);
            elseif InputDataParams.direction == 2;
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_reverse(InputDataParams);
            else
                [Data,Locations,RayParams,BackAzimuths] = build_input_matrices_both_directions(InputDataParams);
            end
            fprintf('done.\n')

            [Data,obj.Locations,obj.RayParams,obj.BackAzimuths] = random_sample(Data,Locations,RayParams,BackAzimuths,0.0);
           

            obj.d=reshape(Data',[],1);
            obj.InputDataParams=InputDataParams;
            
            
        end
        function seis=extractSeis(obj,iseis)
            ind1=(iseis-1)*obj.nptsPerSeis+1;
            ind2=ind1+obj.nptsPerSeis-1;
            seis=obj.d(ind1:ind2);
        end
    end

end
