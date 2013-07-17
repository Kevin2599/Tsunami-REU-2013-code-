function fun = fixit(temp, s1, s2, iterations)% domain, Function, and indexes of region to smooth
for i=1:iterations
   j=s1;
   while j<=s2
       temp(j)=(temp(j-1)+temp(j+1))/2;
       j=j+1;
   end
end
fun=temp;
end
