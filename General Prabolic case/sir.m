n=10000;
y=rand(2,1);
 figure(1)
 hold on
for i=1:n
 x=randi(3);
  if x==1
      y(1)=y(1)/2;
      y(2)=y(2)/2;
       plot(y(1),y(2),'.r','MarkerSize',3)
  end
  
 if x==2
      y(1)=(y(1)+1)/2;
      y(2)=y(2)/2;
       plot(y(1),y(2),'.b','MarkerSize',3)
 end
 
 if x==3
           y(1)=(y(1)+sqrt(1)/2)/2;
      y(2)=(y(2)+sqrt(2)/2)/2;
       plot(y(1),y(2),'.k','MarkerSize',3)
     
 end

end


