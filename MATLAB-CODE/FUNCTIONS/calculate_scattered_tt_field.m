    function [Tscat, Pscat] = calculate_scattered_tt_field(xs,zs)
        
        %xs=-600:10:100;
        %zs=0:1:300;
        
        model.vp=[5.660,7.920,7.280];
        %model.vp=[8.0,8.0,8.0];
        model.hs=[60.0,60.0,180.0];
        
        Tscat=zeros(length(xs),length(zs));
        Pscat=zeros(length(xs),length(zs));
        
        for ix = 1:length(xs);
            for iz = 1:length(zs);
               
                x=xs(ix); z=zs(iz);
                vp=get_v(model,z);
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
        
        function [vp]=get_v(model,z)
            nlay=length(model.vp);
            zref=0;
            for ilay=1:nlay;
                top=zref;
                bot=zref+model.hs(ilay);
                zref=bot;

                if ( z>bot )    %depth below layer
                   %
                elseif (z<top)  %depth above layer
                   %
                else             %depth within layer
                   vp=model.vp(ilay);
                   return
                end

            end
            
            vp=NaN;
            return
            
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
