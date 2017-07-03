function [Tdirect,Tscat,Pscat,xs,zs]=read_in_necessary_files(rp, xs, zs, model)

    rp=rp/111.1;
    
    %xs=(-600:1:100)';
    %zs=0:1:300;
    %[Pscat,xs,zs] = read_table('MIGRA/ILEG_2/tabl
    %[Tscat,~,~] = read_table('MIGRA/ILEG_2/table_time.txt');

    [Tdirect]=calculate_travel_time_field();
    
    [Tscat,Pscat] = calculate_scattered_tt_field();
    
    [XS,~]=meshgrid(xs,zs);
    
    Pscat=-Pscat.*sign(XS');
    
    Pscat=Pscat*111.1;
    
    function [Tscat, Pscat] = calculate_scattered_tt_field()
        
        %xs=-600:10:100;
        %zs=0:1:300;
        
        %model.vp=[5.660,7.920,7.280];
        %model.hs=[60.0,60.0,180.0];
        
        Tscat=zeros(length(xs),length(zs));
        Pscat=zeros(length(xs),length(zs));
        
        for ix = 1:length(xs);
            for iz = 1:length(zs);
               
                x=xs(ix); z=zs(iz);
                [vp,~]=get_v(model,z);
                pmax=1/vp;
                
                dp=0.0001;
                p0=0.0;
                
                    for iter=1:100;
                        %function f(x) = xray - abs(x)

                        [tray,xray]=shoot_ray_up(p0,z,model);
                        [~,xray2]=shoot_ray_up(p0+dp,z,model);
                        dxray=xray2-xray;
                        df_dp=(dxray)/dp;
                        p0=p0 - (xray-abs(x))/df_dp;
                        
                        %if iter>10
                            %fprintf('Warning many iterations needed for x, z, iter = %f %f %d\n',x,z, iter)
                            %fprintf('iter, p0, xray, dxray, dp = %d %f %f %f %f %f\n',iter, p0, xray, dxray, (xray-abs(x))/df_dp, pmax)
                        %end

                        if p0>pmax || imag(tray) ~= 0.0;
                            %if overshoot==1;
                                Tscat(ix,iz)=NaN;
                                Pscat(ix,iz)=NaN;
                                %fprintf('Overshot !\n')
                                TryAlgorithm2=1;
                                break
                            
                        elseif x == 0;
                            Tscat(ix,iz)=tray;
                            Pscat(ix,iz)=p0;
                            TryAlgorithm2=0;
                            break
                        elseif abs((xray-abs(x))/df_dp) < 1.0E-5;
                            tkeep=tray;
                            pkeep=p0;
                            Tscat(ix,iz)=tkeep;
                            Pscat(ix,iz)=pkeep;
                            TryAlgorithm2=0;
                            break
                            
                        elseif isnan(p0)
                            Tscat(ix,iz)=NaN;
                            Pscat(ix,iz)=NaN;
                            TryAlgorithm2=1;
                            break
                        else
                            Tscat(ix,iz)=NaN;
                            Pscat(ix,iz)=NaN;
                            TryAlgorithm2=1;
                        end

                    end
                
                if (TryAlgorithm2==1)
                    ang1=0.0; ang2=90.0; np=100;
                    angles=linspace(ang1,ang2,np);
                    
                    for angIter = 1:4;
                    
                        for angle = angles;
                            p=sind(angle)/vp;
                            if p==pmax;
                                break
                            end
                            [tray,xray]=shoot_ray_up(p,z,model);
                            if xray >= abs(x);
                                %fprintf('x,z,xray,tray = %f, %f, %f, %f \n',x,z,xray,tray)
                                break
                            end
                            tray_last=tray;
                            xray_last=xray;
                            p_last=p;
                        end

                        if p==pmax;
                            Tscat(ix,iz)=NaN;
                            Pscat(ix,iz)=NaN;
                        elseif x == 0;
                            Tscat(ix,iz)=tray;
                            Pscat(ix,iz)=p;
                        else
                            tkeep=interp1([xray_last,xray],[tray_last,tray],abs(x));
                            pkeep=interp1([xray_last,xray],[p_last,p],abs(x));
                            Tscat(ix,iz)=tkeep;
                            Pscat(ix,iz)=pkeep;
                        end
                        
                        if (p==pmax)
                            ang1=angles(end-1);
                            ang2=angles(end);
                            angles=linspace(ang1,ang2,np);
                        else
                            break
                        end
                        
                    end             

                end        
%    
            end
                        
        end
        
        function [T,X]=shoot_ray_up(p,z,model)
            vs=model.vp;
            hs=model.hs;

            nlay=length(hs);
            T=0.0;
            X=0.0;
            zref=0.0;
            for ilay=1:nlay;
                top=zref;
                bot=zref+model.hs(ilay);
                zref=bot;

                if ( z>=bot )    %depth below layer
                    dz=hs(ilay);
                elseif (z<=top)  %depth above layer
                    dz=0.0;
                else             %depth within layer
                    dz=hs(ilay)-(bot-z);    
                end
             
                u=1./vs(ilay);
                
                T=T+u.^2*dz/(u.^2-p.^2).^0.5;
                X=X+p.*dz/(u.^2-p.^2).^0.5;
                
                if X==Inf;
                    error('Error in shoot ray up.')
                end
            end
        end
        
    end

    function [TT]=calculate_travel_time_field()
        
        %model.vs=[3.2,4.4,4.05];
        %model.hs=[60.0,60.0,180.0];

        Tss=zeros(1,length(zs));
        Xss=zeros(1,length(zs));
        i=0;
        for z = zs;
            i=i+1;
            [T,X]=shoot_ray(rp,z,model);

            Tss(i)=T;
            Xss(i)=X;
            Tss(i)=T-rp*X;
            %fprintf('%f %f %f %f\n',z,T,X,Tss(i))
        end
        %clf;
        %plot(Xss,Tss);

        [XS,TSS]=meshgrid(xs,Tss);
        TT = TSS + XS*rp;

        %contourf(xs,-zs,TT); colorbar;
        function [T,X]=shoot_ray(p,z,model)
            vs=model.vs;
            hs=model.hs;

            nlay=length(hs);
            T=0.0;
            X=0.0;
            zref=0.0;
            for ilay=1:nlay;
                top=zref;
                bot=zref+model.hs(ilay);

                if ( z>=bot )    %depth below layer
                    dz=0.0;
                elseif (z<=top)  %depth above layer
                    dz=hs(ilay);
                else             %depth within layer
                    dz=hs(ilay)-(z-top);    
                end

                zref=bot;
                u=1./vs(ilay);
                T=T+u.^2*dz/(u.^2-p.^2).^0.5;
                X=X+p.*dz/(u.^2-p.^2).^0.5;
            end
        end

    end

end
