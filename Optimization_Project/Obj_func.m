function y= Obj_func(x,a,c)

sum=0;
basis=@(xi,ai,ci) ai*(xi/(1-xi/ci));


for i=1:17
    sum=sum+basis(x(i),a(i),c(i));
end

y=sum;

end

