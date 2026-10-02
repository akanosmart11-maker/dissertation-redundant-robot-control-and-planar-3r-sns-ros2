function [q_dot, s_star] = SNS_velocity_p(J,x_dot,Vmin,Vmax,n,m)
    W = eye(n); q_dot_N = zeros(n,1);
    s_star = 0; q_dot_star = zeros(n,1);
    for iter = 1:(n+1)
        JW = J*W;
        if rank(JW,1e-6) < m; break; end
        JWp = pinv(JW);
        q_dot = q_dot_N + JWp*(x_dot - J*q_dot_N);
        if all(q_dot >= Vmin-1e-10) && all(q_dot <= Vmax+1e-10)
            s_star=1; q_dot_star=q_dot; break;
        end
        a=JWp*x_dot; b=q_dot-a;
        sv=ones(n,1);
        for i=1:n
            if abs(a(i))<1e-10
                if b(i)<Vmin(i)-1e-10||b(i)>Vmax(i)+1e-10
                    sv(i)=0;
                end
            else
                sL=(Vmin(i)-b(i))/a(i); sU=(Vmax(i)-b(i))/a(i);
                if sL>sU; tmp=sL; sL=sU; sU=tmp; end
                if sU<0||sL>1
                    sv(i)=0;
                else
                    sv(i)=min(1,max(0,sU));
                end
            end
        end
        [sk,jc]=min(sv);
        if sk>s_star
            s_star=sk;
            q_dot_star=max(Vmin,min(Vmax,...
                q_dot_N+JWp*(sk*x_dot-J*q_dot_N)));
        end
        W(jc,jc)=0;
        if q_dot(jc)>Vmax(jc)
            q_dot_N(jc)=Vmax(jc);
        else
            q_dot_N(jc)=Vmin(jc);
        end
    end
    q_dot = max(Vmin,min(Vmax,q_dot_star));
end