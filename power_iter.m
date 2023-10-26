function max_egival=power_iter(A,At,image)
max_iter=100;
b0=rand(size(image));
b0=b0./normNdMatrix(b0,2);
iter=0;
cond=1;
while(cond)
    iter=iter+1;
    if(iter==1)
        b_old=b0;
    end
    b_new=At(A(b_old));  
    b_old_temp=b_old(:);
    u_new=b_old_temp.'*b_new(:);
    b_new=b_new./normNdMatrix(b_new,2);    
    if(iter>1)
        if((u_new-u_old)/u_new<1e-3)||(iter>=max_iter)
            cond=0;
        end   
    end
    b_old=b_new;
    u_old=u_new;
end
max_egival=u_new;
iter

function norm_val = normNdMatrix(x,n)
norm_val_temp = x.^n;
norm_val=sum(norm_val_temp(:));
norm_val=norm_val^(1/n);
