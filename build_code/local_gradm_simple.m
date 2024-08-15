function [S,A,adh,adL]=local_gradm_simple(L)
% function  [S,A]=local_gradm(X,Y,L)
%   computes slope and aspect for surface
%   looks at all eight neightbors and finds steepest descent
%
% not using preexisting products
%   "gradientm" makes too many assumption so gives wrong result
%   "gradient" only computes along dimentions, not all 8 neighbors

% get size of grids
[m,n]=size(L);

% create X, Y grids with binsize=1
%[X,Y]=meshgrid(1:m,1:n);

% create blank matricies for slope and aspect
S=zeros(size(L));
A=zeros(size(L));
adh=zeros(size(L));
adL=zeros(size(L));
%
% set up x & y bin distances (assumes all bins 1x1)
        % get x-dir bin distances
        dXset(1)=1; % X(i,j+1)-X(i,j);  
        dXset(2)=1; % X(i+1,j+1)-X(i,j);
        dXset(3)=0; % X(i+1,j)-X(i,j);
        dXset(4)=-1; % X(i+1,j-1)-X(i,j);
        dXset(5)=-1; % X(i,j-1)=X(i,j);
        dXset(6)=-1; % X(i-1,j-1)-X(i,j);
        dXset(7)=0; % X(i-1,j)=X(i,j);
        dXset(8)=1; % X(i-1,j+1)X(i,j);
        % get y-dir bin distances
        dYset(1)=0; % Y(i,j+1-Y(i,j)-;  
        dYset(2)=1; % Y(i+1,j+1)-Y(i,j);
        dYset(3)=1; % Y(i+1,j)-Y(i,j)-;
        dYset(4)=1; % Y(i+1,j-1)-Y(i,j);
        dYset(5)=0; % Y(i,j-1)-Y(i,j);
        dYset(6)=-1; % Y(i-1,j-1)-Y(i,j);
        dYset(7)=-1; % Y(i-1,j)-Y(i,j);
        dYset(8)=-1; % Y(i-1,j+1)-Y(i,j);
        % then precalculate bin-to-bin distances
        dHset=sqrt(dXset.^2 + dYset.^2);
        % also pre-calc bin-to-bin directions
        Aset=atan2d(dYset,dXset);
%
% step though points
% interior points consider all eight neighbors
for i=2:m-1
    for j=2:n-1
        % get elevation differences
        dLset(1)=L(i,j)-L(i,j+1);  
        dLset(2)=L(i,j)-L(i+1,j+1);
        dLset(3)=L(i,j)-L(i+1,j);
        dLset(4)=L(i,j)-L(i+1,j-1);
        dLset(5)=L(i,j)-L(i,j-1);
        dLset(6)=L(i,j)-L(i-1,j-1);
        dLset(7)=L(i,j)-L(i-1,j);
        dLset(8)=L(i,j)-L(i-1,j+1);

        % find the max elevation change and its direction
        dLdHset=dLset./dHset;
        dLdH=max(dLdHset);
        dirset=find(dLdHset==dLdH);
        k=length(dirset); 
        kcon=max(diff(dirset)); % changed from k to dirset 8/22/22 kgb
        if k>2
            display(dirset)
            fprintf('dLdH = %f; k %d; kcon %d \n',dLdH,k,kcon)
        end
        %if dLdH<0.001
        %    % if slope very small, make direction totally random
        %    v=randi(8,1);
        %    dir=v;
        %    % otherwise look for steepest line of descent
        %elseif k==1
        if k==1
            dir=dirset;
        elseif kcon<=1
            % this method better if small k with contiguous dirset values
            % and applies only when slope large enough to be meaningful
            % (slope criteria added 8/22/22 kgb)
            if dLdH>0.0001
                v=round(mean(dirset));
                dir=v;
            else
                v=randi(k,1);
                dir=dirset(v);
            end
        else
        	% should use this method if k>>3 or dirset values not
        	% contiguous
        	v=randi(k,1);
        	dir=dirset(v);
        end
        
        dh=dHset(dir);
        
        adh(i,j)=dh;
        adL(i,j)=dLset(dir);
        S(i,j)=dLdH;
        A(i,j)=Aset(dir);
        %display(dLset)
    end
end
% edge points consider limited neighbors
% first row, all columns
for i=1
    for j=1:n
        dLset=NaN*ones(8,1)';
        if j<n
            dLset(1)=L(i,j)-L(i,j+1);  
            dLset(2)=L(i,j)-L(i+1,j+1);
        end
        dLset(3)=L(i,j)-L(i+1,j);
        if j>1
        dLset(4)=L(i,j)-L(i+1,j-1);
        dLset(5)=L(i,j)-L(i,j-1);
        end
        dLdHset=dLset./dHset;
	% changing from 
	%dLdH=max(dLdHset,[],"omitnan");
	% to
	dLdH=max(dLdHset(~isnan(dLdHset)));
	%
        dirset=find(dLdHset==dLdH);
        k=length(dirset);
        kcon=max(diff(dirset));
        %if dLdH<0.001
        %    % if slope very small, make direction totally random
        %    v=randi(8,1);
        %    dir=v;
        %    % otherwise look for steepest line of descent
        %elseif k==1
        if k==1
            dir=dirset;
        elseif kcon<=1
            % this method better if small k with contiguous dirset values
            v=round(mean(dirset));
            dir=v;
        else
        	% should use this method if k>>3 or dirset values not
        	% contiguous
        	v=randi(k,1);
        	dir=dirset(v);
        end
        %display(dir)
        dh=dHset(dir);
        
        adh(i,j)=dh;
        adL(i,j)=dLset(dir);
        S(i,j)=dLdH;
        A(i,j)=Aset(dir);
        %display(dLset)
    end
end
% last row, all columns
for i=n
    for j=1:n
        dLset=NaN*ones(8,1)';
        if j<n
            dLset(1)=L(i,j)-L(i,j+1);  
            dLset(8)=L(i,j)-L(i-1,j+1);
        end
        dLset(7)=L(i,j)-L(i-1,j);
        if j>1
            dLset(5)=L(i,j)-L(i,j-1);
            dLset(6)=L(i,j)-L(i-1,j-1);            
        end
        dLdHset=dLset./dHset;
        % changing from
	% dLdH=max(dLdHset,[],"omitnan");
	% to 
	dLdH=max(dLdHset(~isnan(dLdHset)));
	%
        dirset=find(dLdHset==dLdH);
        k=length(dirset);
        kcon=max(diff(dirset));
        %if dLdH<0.001
        %   % if slope very small, make direction totally random
        %    v=randi(8,1);
        %    dir=v;
        %    % otherwise look for steepest line of descent
        %elseif k==1
        if k==1
            dir=dirset;
        elseif kcon<=1
            % this method better if small k with contiguous dirset values
            v=round(mean(dirset));
            dir=v;
        else
        	% should use this method if k>>3 or dirset values not
        	% contiguous
        	v=randi(k,1);
        	dir=dirset(v);
        end
      
        dh=dHset(dir);
        
        adh(i,j)=dh;
        adL(i,j)=dLset(dir);
        S(i,j)=dLdH;
        A(i,j)=Aset(dir);
        
        %display(dLset)
    end
end
% more edges: first column, inner rows
for i=2:m-1
    for j=1
        dLset=NaN*ones(8,1)';
        % get elevation differences
        dLset(1)=L(i,j)-L(i,j+1);  
        dLset(2)=L(i,j)-L(i+1,j+1);
        dLset(3)=L(i,j)-L(i+1,j);
        dLset(7)=L(i,j)-L(i-1,j);
        dLset(8)=L(i,j)-L(i-1,j+1);

        % find the max elevation change and its direction
        dLdHset=dLset./dHset;
        dLdH=max(dLdHset);
        dirset=find(dLdHset==dLdH);
        k=length(dirset); 
        kcon=max(diff(dirset));
        %if dLdH<0.001
        %    % if slope very small, make direction totally random
        %    v=randi(8,1);
        %    dir=v;
        %    % otherwise look for steepest line of descent
        %elseif k==1
        if k==1
            dir=dirset;
        elseif kcon<=1
            % this method better if small k with contiguous dirset values
            v=round(mean(dirset));
            dir=v;
        else
        	% should use this method if k>>3 or dirset values not
        	% contiguous
        	v=randi(k,1);
        	dir=dirset(v);
        end
        
        dh=dHset(dir);
        
        adh(i,j)=dh;
        adL(i,j)=dLset(dir);
        S(i,j)=dLdH;
        A(i,j)=Aset(dir);
        %display(dLset)
    end
end

% more edges: last column, inner rows
for i=2:m-1
    for j=n
        dLset=NaN*ones(8,1)';
        % get elevation differences
        dLset(3)=L(i,j)-L(i+1,j);
        dLset(4)=L(i,j)-L(i+1,j-1);
        dLset(5)=L(i,j)-L(i,j-1);
        dLset(6)=L(i,j)-L(i-1,j-1);
        dLset(7)=L(i,j)-L(i-1,j);

        % find the max elevation change and its direction
        dLdHset=dLset./dHset;
        dLdH=max(dLdHset);
        dirset=find(dLdHset==dLdH);
        k=length(dirset); 
        kcon=max(diff(dirset));
        %if dLdH<0.001
        %    % if slope very small, make direction totally random
        %    v=randi(8,1);
        %    dir=v;
        %    % otherwise look for steepest line of descent
        %elseif k==1
        if k==1
            dir=dirset;
        elseif kcon<=1
            % this method better if small k with contiguous dirset values
            v=round(mean(dirset));
            dir=v;
        else
        	% should use this method if k>>3 or dirset values not
        	% contiguous
        	v=randi(k,1);
        	dir=dirset(v);
        end
        
        dh=dHset(dir);
        
        adh(i,j)=dh;
        adL(i,j)=dLset(dir);
        S(i,j)=dLdH;
        A(i,j)=Aset(dir);
        %display(dLset)
    end
end



