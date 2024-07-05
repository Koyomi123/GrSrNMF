function [U,V,Z,W,D,Zsnmf,obj,Q]=DESNMF(W_Cube,GT_Cube,U0,V0,alpha,beta,gama,type)

for t = 1:length(W_Cube)
    U{t} = U0{t};
    V{t} = V0{t};
end
[W,D]=computeSimilarity(W_Cube,type);
[~,Zsnmf(:,1)]=max(V{1}');
final_community_assignment = Zsnmf(:,1);
Q{1} = modularity(W{1},final_community_assignment);
for t = 2:length(W_Cube)
    Z{t-1}=rand(size(GT_Cube{t-1},2),size(GT_Cube{t},2));
    i = 0;
    Jnew = 10;
    Jold = 1;
%     U{t} = U0{t};
%     V{t} = V0{t};
    while abs(Jnew - Jold)/Jold > 1e-6
        i = i + 1;
%         for s=1:t-1
%             LL{t}= theta^(t-s)*L{s};
%             WW{t}= theta^(t-s)*W{s};
%             D{t}=LL{t}+WW{t};
%         end
        U{t}=U{t}.*(((2*gama*V{t}*U{t}'*V{t})+(W_Cube{t}*V{t}))./max(realmin,((U{t}*V{t}'*V{t})+(2*gama*U{t}*V{t}'*V{t}))));
        V{t}=V{t}.*((2*gama*U{t}*V{t}'*U{t}+W_Cube{t}'*U{t}+alpha*V{t-1}*Z{t-1}+beta*W{t-1}*V{t})./max(realmin,(V{t}*U{t}'*U{t}+alpha*V{t}+beta*D{t-1}*V{t}+2*gama*V{t}*U{t}'*U{t})));
        Z{t-1} = Z{t-1}.*((V{t-1}'*V{t})./max(realmin,(V{t-1}'*V{t-1}*Z{t-1})));
        Jold = Jnew;
        L1=norm(W_Cube{t} - U{t}*V{t}','fro')^2;
        L2=alpha*norm(V{t-1}*Z{t-1}-V{t},'fro')^2;
        L3=gama*norm(U{t}*V{t}'-V{t}*U{t}','fro')^2;%ÐÂµÄÏî
        L4=beta*trace(V{t}'*(D{t-1}-W{t-1})*V{t});

        Jnew = L1+L2+L3+L4;
            obj(i,t)=Jnew;
            if i>800
                break;
            end
    end
     [~,Zsnmf(:,t)]=max(V{t}');
     final_community_assignment = Zsnmf(:,t);
     Q{t} = modularity(W{t},final_community_assignment);
     
end