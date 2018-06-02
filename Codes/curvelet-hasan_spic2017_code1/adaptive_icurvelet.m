function x = adaptive_icurvelet(C, is_real, N1, N2,scale_locations, nbangles)
% Inputs
%   C           Cell array containing curvelet coefficients (see
%               description in fdct_wrapping.m)
%   is_real     Type of the transform
%                   0: complex-valued curvelets
%                   1: real-valued curvelets
%               [default set to 0]
%   M, N        Size of the image to be recovered (not necessary if finest
%               = 2)
%
% Outputs
%   x           M-by-N matrix
%
% See also fdct_wrapping.m
%
% By Hasan Al-Marzouqi as a modification of the wrapping curvelet code.


%Default is 4 levels:
if nargin<5, scale_locations=floor([(((4/3*N1-4/6*N1)/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)  ((4/3*N1-4/6*N1)/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)    (4/3*N1-4/6*N1)/2+(4/6*N1)  4/3*N1
                                (((4/3*N2-4/6*N2)/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)  ((4/3*N2-4/6*N2)/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)    (4/3*N2-4/6*N2)/2+(4/6*N2)  4/3*N2]);end;

if nargin<6,nbangles_coarse=16; const=4; nbangles=[1, nbangles_coarse .* 2.^(ceil((const-(const:-1:2))/2))];nbangles=[nbangles;nbangles];end;




scale_locations=ceil(scale_locations);

numscales=size(scale_locations,2);



M1=scale_locations(1,end);
M2=scale_locations(2,end);

M1S=((M1-1)/4);
M2S=((M2-1)/4);

    
bigN1 = 2*floor(2*M1S)+1;
bigN2 = 2*floor(2*M2S)+1;
X = zeros(bigN1,bigN2);
%finding size of smoothigwindow:
smoothwinsize=zeros(2,numscales);
smoothwinsize(:,end)=(scale_locations(:,end)-[N1;N2]);
smoothwinsize(:,end)=(scale_locations(:,end)-[N1;N2]);

if mod(smoothwinsize(1,end),2)==0
    smoothwinsize(1,end)=smoothwinsize(1,end)-1;
end

if mod(smoothwinsize(2,end),2)==0
    smoothwinsize(2,end)=smoothwinsize(2,end)-1;               
end

    %changing scale_locations from original coordinates to how far each level is
%from the origin:
scale_locations(1,1:end)=scale_locations(1,1:end)-floor(bigN1/2)*ones(1,size(scale_locations,2));
scale_locations(2,1:end)=scale_locations(2,1:end)-floor(bigN2/2)*ones(1,size(scale_locations,2));

for i=numscales-1:-1:1
    delta_ip=scale_locations(:,i+1)-scale_locations(:,i);
    if i~=1
        delta_in=scale_locations(:,i)-scale_locations(:,i-1);
    else delta_in=scale_locations(:,1)-1;
    end 
    delta_i0=scale_locations(:,i)-1;
    Thsmooth=1.25; 
    smoothwinsize(:,i)=min(max(delta_ip/Thsmooth,delta_in/Thsmooth),min(delta_i0,delta_ip));                           
end








M1S=floor(smoothwinsize(1,numscales));
M2S=floor(smoothwinsize(2,numscales));




window_length_1 = floor(2*M1S)-floor(M1S)  -1; 
window_length_2 = floor(2*M2S)-floor(M2S) -1;



coord_1 = 0:(1/window_length_1):1;
coord_2 = 0:(1/window_length_2):1;    
[wl_1,wr_1] = fdct_wrapping_window(coord_1);
[wl_2,wr_2] = fdct_wrapping_window(coord_2);    
lowpass_1 = [wl_1, ones(1,bigN1-length(wl_1)-length(wr_1)), wr_1];
lowpass_2 = [wl_2, ones(1,bigN2-length(wl_2)-length(wr_2)), wr_2];
lowpass = lowpass_1'*lowpass_2;  
scales = numscales:-1:2;   
    

% Loop: pyramidal reconstruction
startindex1=0;
startindex2=0;
Xj_topleft_1 = 1;
Xj_topleft_2 = 1;
Xj=X;
Xjnew=X;
for j = scales,
      %j;
   sizecurrentlevel=size(Xjnew);
   newN1=size(Xj,1);
   newN2=size(Xj,2);
   M1=scale_locations(1,j-1);
   M2=scale_locations(2,j-1);


   Xlow_index_1 = ((-floor(M1+1/2*smoothwinsize(1,j-1))):floor(M1+1/2*smoothwinsize(1,j-1))) + floor(newN1/2) + 1;
   Xlow_index_2 = ((-floor(M2+1/2*smoothwinsize(2,j-1))):floor(M2+1/2*smoothwinsize(2,j-1))) + floor(newN2/2) + 1;


   Xjnew=zeros(length(Xlow_index_1),length(Xlow_index_2));
   newN1=size(Xjnew,1);
   newN2=size(Xjnew,2);
    
      M1S=floor(smoothwinsize(1,j-1));
      M2S=floor(smoothwinsize(2,j-1));


    window_length_1 = floor(2*M1S) - floor(M1S) - 1;
    window_length_2 = floor(2*M2S) - floor(M2S) - 1;
    coord_1 = 0:(1/window_length_1):1;
    coord_2 = 0:(1/window_length_2):1;
    [wl_1,wr_1] = fdct_wrapping_window(coord_1);
    [wl_2,wr_2] = fdct_wrapping_window(coord_2);
    lowpass_1 = [wl_1, ones(1,newN1-length(wr_1)-length(wr_1)), wr_1];
    lowpass_2 = [wl_2, ones(1,newN2-length(wl_2)-length(wr_2)), wr_2];      
    lowpass_next = lowpass_1'*lowpass_2;
    hipass = sqrt(1 - lowpass_next.^2);   
    % Loop: angles
    l = 0;
    nbquadrants = 2 + 2*(~is_real);
    
    if nbangles(1,j)~=nbangles(2,j)
        notequalangular=1;
    else     notequalangular=0;
    end

    vert=floor((sizecurrentlevel(1))/2);
    horiz=floor((sizecurrentlevel(2))/2);
    
    for quadrant = 1:nbquadrants
            if nbangles(1,j)<nbangles(2,j)
                    tempquadrant=2;
                     myfourMhoriz= horiz * (mod(quadrant,2)==1) + vert * (mod(quadrant,2)==0);        
                     myfourMvert= vert * (mod(quadrant,2)==1) + horiz * (mod(quadrant,2)==0);        
                     M_horiz=(myfourMhoriz/4);
                     M_vert=(myfourMvert/4);
                    nbangles_perquad = nbangles(tempquadrant,j);
                    if mod(nbangles_perquad,2),            
                        wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
                        wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
                         wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
                    else
                        wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
                        wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
                        wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];                        
                    end;
                    wedge_endpoints = wedge_ticks(2:2:(end-1)); % integers
            
                    first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);      

                    C2left = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1)+ 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));

                    C1left = C2left / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
                    
                    C2right = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
                    
                    C1right = -C2right * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);

                    
            else
                    tempquadrant=1;                    
                    myfourMhoriz= horiz * (mod(quadrant,2)==1) + vert * (mod(quadrant,2)==0);        
                     myfourMvert= vert * (mod(quadrant,2)==1) + horiz * (mod(quadrant,2)==0);        
                     M_horiz=(myfourMhoriz/4);
                     M_vert=(myfourMvert/4);

                    nbangles_perquad = nbangles(tempquadrant,j);
                    if mod(nbangles_perquad,2),            
                        wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
                        wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
                         wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
                    else
                        wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
                        wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
                        wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];
                        
                    end;
                    wedge_endpoints = wedge_ticks(2:2:(end-1)); % integers
            
                    first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles_perquad) + 1);      

                    C2left = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1)+ 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));

                    C1left = C2left / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
                    
                    C2right = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
                    
                    C1right = -C2right * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);


                    
            end                                                 
                                    
         nbangles_perquad = nbangles(quadrant,j);
% 
%         
          myfourMhoriz= horiz * (mod(quadrant,2)==1) + vert * (mod(quadrant,2)==0);        
          myfourMvert= vert * (mod(quadrant,2)==1) + horiz * (mod(quadrant,2)==0);        
          M_horiz=(myfourMhoriz/4);
          M_vert=(myfourMvert/4);




        if mod(nbangles_perquad,2),
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right(end:-1:1)];
        else
            wedge_ticks_left = round((0:(1/(2*nbangles_perquad)):.5)*2*floor(4*M_horiz) + 1);
            wedge_ticks_right = 2*floor(4*M_horiz) + 2 - wedge_ticks_left;
            wedge_ticks = [wedge_ticks_left, wedge_ticks_right((end-1):-1:1)];
        end;
        wedge_endpoints = wedge_ticks(2:2:(end-1));         % integers
        wedge_midpoints = (wedge_endpoints(1:(end-1)) + wedge_endpoints(2:end))/2;

        % Left corner wedge        
        l = l+1;
      if quadrant==1
          otherq=2;
      else otherq=1;
      end

         first_wedge_endpoint_vert = round(2*floor(4*M_vert)/(2*nbangles(otherq,j))+1 );
        
        

         length_corner_wedge=floor(4*M_vert)-floor(0.5*smoothwinsize(quadrant,j-1));
        
        Y_corner = 1:length_corner_wedge;
        [XX,YY] = meshgrid(1:(2*floor(4*M_horiz)+1),Y_corner);
        width_wedge = wedge_endpoints(2) + wedge_endpoints(1) - 1;
        slope_wedge = (floor(4*M_horiz) + 1 - wedge_endpoints(1))/floor(4*M_vert);
        left_line = round(2 - wedge_endpoints(1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;

        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);
                                                            % not integers
                                                            % in general
        coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(2) - wedge_endpoints(1)) * ...
            (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) + 1;
        
        if notequalangular==1
            coord_corner = C1left + C2left * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        else
            C2 = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1) + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
            C1 = C2 / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);
            coord_corner = C1 + C2 * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        end


        wl_left = fdct_wrapping_window(coord_corner);
        [~,wr_right] = fdct_wrapping_window(coord_right);

        
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(numel(C{j}{l}));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1 
            x = C{j}{l} + sqrt(-1)*C{j}{l+sum(nbangles(1,j)+nbangles(2,j))};                
            wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(numel(x))/sqrt(2);
            wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;

       
        wrapped_data = wrapped_data .* (wl_left .* wr_right);

        % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            Xj(row,admissible_cols) = Xj(row,admissible_cols) + wrapped_data(new_row,:);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
        end;
                                                        
        
        
        
        
        
        % Regular wedges

        length_wedge = floor(4*M_vert)-floor(0.5*smoothwinsize(quadrant,j-1));


        Y = 1:length_wedge;
        first_row = floor(4*M_vert)+2-ceil((length_wedge+1)/2)+...
            mod(length_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        for subl = 2:(nbangles_perquad-1);
            l = l+1;
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert);
            left_line = round(wedge_endpoints(subl-1) + slope_wedge*(Y - 1));
            [wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                wrapped_XX(new_row,:) = XX(row,cols);
                wrapped_YY(new_row,:) = YY(row,cols);
            end;
            slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(subl-1))/floor(4*M_vert);
            mid_line_left = wedge_midpoints(subl-1) + slope_wedge_left*(wrapped_YY - 1);
            coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl) - wedge_endpoints(subl-1)) * ...
                (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
            slope_wedge_right = ((floor(4*M_horiz)+1) - wedge_midpoints(subl))/floor(4*M_vert);
            mid_line_right = wedge_midpoints(subl) + slope_wedge_right*(wrapped_YY - 1);
            coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(subl+1) - wedge_endpoints(subl)) * ...
                (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);
            wl_left = fdct_wrapping_window(coord_left);
            [~,wr_right] = fdct_wrapping_window(coord_right);
            
            

            switch is_real
             case 0
              wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(numel(C{j}{l}));
              wrapped_data = rot90(wrapped_data,(quadrant-1));
             case 1

                   x = C{j}{l} + sqrt(-1)*C{j}{l+sum(nbangles(1,j)+nbangles(2,j))};
                   
                  wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(numel(x))/sqrt(2);
                  wrapped_data = rot90(wrapped_data,(quadrant-1));
            end;
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            
            % Unwrapping data
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
                new_row = 1 + mod(row - first_row, length_wedge);
                Xj(row,cols) = Xj(row,cols) + wrapped_data(new_row,:);
            end;

        end;    % for subl
        
        % Right corner wedge
        l = l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        [wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
        
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);        
        end;
        
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert)+1 - wrapped_YY);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY-1)/floor(4*M_vert)) - 1;

        if notequalangular==1
            coord_corner = C1right + C2right * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));

        else
            C2 = -1/(2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1 + 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
            C1 = -C2 * (2*(floor(4*M_horiz))/(wedge_endpoints(end) - 1) - 1);
            coord_corner = C1 + C2 * (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert)))) ./ ...
            ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert)));
            
            
        end
        
        wl_left = fdct_wrapping_window(coord_left);
        
        [~,wr_right] = fdct_wrapping_window(coord_corner);
        
        switch is_real
         case 0
          wrapped_data = fftshift(fft2(ifftshift(C{j}{l})))/sqrt(numel(C{j}{l}));
          wrapped_data = rot90(wrapped_data,(quadrant-1));
         case 1
               x = C{j}{l} + sqrt(-1)*C{j}{l+sum(nbangles(1,j)+nbangles(2,j))};
             
             
             
          wrapped_data = fftshift(fft2(ifftshift(x)))/sqrt(numel(x))/sqrt(2);
          wrapped_data = rot90(wrapped_data,(quadrant-1));
        end;
        
        wrapped_data = wrapped_data .* (wl_left .* wr_right);                

         % Unwrapping data
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);                                    
            Xj(row,fliplr(admissible_cols)) = Xj(row,fliplr(admissible_cols)) + wrapped_data(new_row,end:-1:1);
                                % We use the following property: in an assignment
                                % A(B) = C where B and C are vectors, if
                                % some value x repeats in B, then the
                                % last occurrence of x is the one
                                % corresponding to the eventual assignment.
        end;
                                                                           
        Xj = rot90(Xj);

    end;    % for quadrant

       [a,b]=size(Xj);
       Xj = Xj .* lowpass;       
      
    
       M1=scale_locations(1,j-1);       
       M2=scale_locations(2,j-1);   
%        
       
       Xj_index1 = ((-floor(M1+1/2*smoothwinsize(1,j-1))):floor(M1+1/2*smoothwinsize(1,j-1))) + floor(a/2) + 1;
       Xj_index2 = ((-floor(M2+1/2*smoothwinsize(2,j-1))):floor(M2+1/2*smoothwinsize(2,j-1))) + floor(b/2) + 1;


       
       Xj(Xj_index1, Xj_index2) = Xj(Xj_index1, Xj_index2) .* hipass;

    if j==scales(1)
        loc_1=Xj_topleft_1 + (1:a)-1;
        loc_2=Xj_topleft_2 + (1:b)-1;
    else
        loc_1=startindex1+(min(storedXj_index1):max(storedXj_index1));
        loc_2=startindex2+(min(storedXj_index2):max(storedXj_index2));        
        startindex1=min(loc_1)-1;
        startindex2=min(loc_2)-1;
    end
    
   
    storedXj_index1=Xj_index1;
    storedXj_index2=Xj_index2;

       
    X(loc_1,loc_2) = X(loc_1,loc_2) + Xj;
               
       
    lowpass = lowpass_next;  
     
    Xj=Xjnew;       
end;    % for j


if is_real
    Y = X;
    X = rot90(X,2);
    X = X + conj(Y);
end




Xj = fftshift(fft2(ifftshift(C{1}{1})))/sqrt(numel(C{1}{1}));


loc_1=startindex1+(min(storedXj_index1):max(storedXj_index1));
loc_2=startindex2+(min(storedXj_index2):max(storedXj_index2));        

X(loc_1,loc_2) = X(loc_1,loc_2) + Xj .* lowpass;





% Finest level
    % Folding back onto N1-by-N2 matrix
shift_1=floor((bigN1-N1)/2);
shift_2=floor((bigN2-N2)/2); 
Y = X(:,(1:N2)+shift_2);    
Y(:,N2-shift_2+(1:shift_2)) = Y(:,N2-shift_2+(1:shift_2)) + X(:,1:shift_2);    
Y(:,1:shift_2) = Y(:,1:shift_2) + X(:,N2+shift_2+(1:shift_2));
X = Y((1:N1)+shift_1,:);
X(N1-shift_1+(1:shift_1),:) = X(N1-shift_1+(1:shift_1),:) + Y(1:shift_1,:);
X(1:shift_1,:) = X(1:shift_1,:) + Y(N1+shift_1+(1:shift_1),:);



x = fftshift(ifft2(ifftshift(X)))*sqrt(numel(X));

if is_real, x = real(x); end;


end

function [wl,wr] = fdct_wrapping_window(x)

% fdct_wrapping_window.m - Creates the two halves of a C^inf compactly supported window
%
% Inputs
%   x       vector or matrix of abscissae, the relevant ones from 0 to 1
%
% Outputs
%   wl,wr   vector or matrix containing samples of the left, resp. right
%           half of the window
%
% Used at least in fdct_wrapping.m and ifdct_wrapping.m
%
% By Laurent Demanet, 2004

wr = zeros(size(x));
wl = zeros(size(x));
x(abs(x) < 2^-52) = 0;
wr((x > 0) & (x < 1)) = exp(1-1./(1-exp(1-1./x((x > 0) & (x < 1)))));
wr(x <= 0) = 1;
wl((x > 0) & (x < 1)) = exp(1-1./(1-exp(1-1./(1-x((x > 0) & (x < 1))))));
wl(x >= 1) = 1;
normalization = sqrt(wl.^2 + wr.^2);
wr = wr ./ normalization;
wl = wl ./ normalization;
end