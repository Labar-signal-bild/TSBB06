function image_rep( im, m ,N)
n= m-5;


A = im2col(im,[N N],'distinct');
size(A)

C = A*A';        
[e l] = eig(C);               
[PM p] = sort(diag(l),'descend'); 
PC = e(:,p);                      

[PC S] = svd(A);                  
PM = diag(S);  

figure;
subplot(2,1,1);plot(PM,'o');
subplot(2,1,2);plot(log(PM),'o');

figure;colormap('gray');
for M=n:m,
    c=PC(:,1:M)'*A;    %Compute coordinates from blocks
    Arec=PC(:,1:M)*c;  %Reconstruct blocks from coordinates
    imrec=col2im(Arec,[N N],size(im),'distinct');  %Reshape into image
    subplot(2,3,M-n+1);imagesc(imrec);axis('off');     %Display image
    title(sprintf('%d principal components',M));   %Set title
end


end

