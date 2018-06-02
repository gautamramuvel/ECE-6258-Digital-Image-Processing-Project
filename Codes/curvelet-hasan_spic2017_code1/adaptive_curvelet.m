function C = adaptive_curvelet(x, is_real,scale_locations, nbangles)
%Function to take the adaptive curvelet tiling for a given image.
% C = fdct_wrapping3_smooth(x, is_real, finest, scale_locations, nbangles,smooth,Thsmooth)
% x: input image
% is_real= 1 (Using only two quadrant).
% finest=1 (Using outer priodic extension window).
% scale_locations: 2xJ vector of horizontal and vertical scale locations.
% nbangles: vector of angular devisions per scale/quadrant pair.
% Example:
% scale_locations =[384   426   512   682;384   426   512   682];
% nbangles = [1    16    32    32;1    16    32    32];
% C=adaptive_curvelet(img,1,1,scale_locations,nbangles);
%
% Written by Hasan Al-Marzouqi, as a modification of the curvelet wrapping
% code (fdct_wrapping).


X = fftshift(fft2(ifftshift(x)))/sqrt(numel(x));
[N1,N2] = size(X);


%testing for errors in scale_locations:
if nargin < 2, is_real = 1; end;

%Default is four levels:
if nargin<3, scale_locations=floor([(((4/3*N1-4/6*N1)/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)  ((4/3*N1-4/6*N1)/2+(4/6*N1)-(4/6*N1))/2+(4/6*N1)    (4/3*N1-4/6*N1)/2+(4/6*N1)  4/3*N1
                                (((4/3*N2-4/6*N2)/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)  ((4/3*N2-4/6*N2)/2+(4/6*N2)-(4/6*N2))/2+(4/6*N2)    (4/3*N2-4/6*N2)/2+(4/6*N2)  4/3*N2]);end
if nargin<4, nbangles_coarse=16; const=4; nbangles=[1, nbangles_coarse .* 2.^(ceil((const-(const:-1:2))/2))];nbangles=[nbangles;nbangles];end;



scale_locations=ceil(scale_locations);

numscales=size(scale_locations,2);

C = cell(1,numscales);

C{1}=cell(1,1);

for j = 2:numscales
  C{j} = cell(1,sum(nbangles(:,j)));
end;

M1=scale_locations(1,end);
M2=scale_locations(2,end);

M1S=((M1-1)/4);
M2S=((M2-1)/4);

    % Initialization: smooth periodic extension of high frequencies
bigN1 = 2*floor(2*M1S)+1;
bigN2 = 2*floor(2*M2S)+1;
  




equiv_index_1 = 1+mod(floor(N1/2)-floor(2*M1S)+(1:bigN1)-1,N1);
equiv_index_2 = 1+mod(floor(N2/2)-floor(2*M2S)+(1:bigN2)-1,N2);        
    
X = X(equiv_index_1,equiv_index_2);        


    

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
Thsmooth=1.25;
for i=numscales-1:-1:1
         delta_ip=scale_locations(:,i+1)-scale_locations(:,i);
        if i~=1
            delta_in=scale_locations(:,i)-scale_locations(:,i-1);
        else delta_in=scale_locations(:,1)-1;
        end 
        delta_i0=scale_locations(:,i)-1;
        smoothwinsize(:,i)=min(max(delta_ip/Thsmooth,delta_in/Thsmooth),min(delta_i0,delta_ip));                            
end



if is_real==0
    nbangles=[nbangles;nbangles];
    smoothwinsize=[smoothwinsize;smoothwinsize];
end
        
M1S=floor(smoothwinsize(1,numscales));
M2S=floor(smoothwinsize(2,numscales));


window_length_1 = floor(2*M1S)-floor(M1S) - 1;
window_length_2 = floor(2*M2S)-floor(M2S) - 1;

coord_1 = 0:(1/window_length_1):1;
coord_2 = 0:(1/window_length_2):1;

[wl_1,wr_1] = fdct_wrapping_window(coord_1);    
[wl_2,wr_2] = fdct_wrapping_window(coord_2);    
lowpass_1 = [wl_1, ones(1,bigN1-length(wl_1)-length(wr_1)), wr_1];
lowpass_2 = [wl_2, ones(1,bigN2-length(wl_2)-length(wr_2)), wr_2];        
lowpass = lowpass_1'*lowpass_2;
Xlow = X .* lowpass;    
scales = numscales:-1:2;    





for j = scales,       
       newN1=size(Xlow,1);
       newN2=size(Xlow,2);
       M1=scale_locations(1,j-1);
       M2=scale_locations(2,j-1);
       Xhi = Xlow;           

       Xlow_index_1 = ((-floor(M1+1/2*smoothwinsize(1,j-1))):floor(M1+1/2*smoothwinsize(1,j-1))) + floor(newN1/2) + 1;
       Xlow_index_2 = ((-floor(M2+1/2*smoothwinsize(2,j-1))):floor(M2+1/2*smoothwinsize(2,j-1))) + floor(newN2/2) + 1;
       Xlow = Xlow(Xlow_index_1, Xlow_index_2);       
     
      
       newN1=size(Xlow,1);
       newN2=size(Xlow,2);
       
       M1S=floor(smoothwinsize(1,j-1));
       M2S=floor(smoothwinsize(2,j-1));
                     
       window_length_1 = floor(2*M1S) - floor(M1S) - 1;
       window_length_2 = floor(2*M2S) - floor(M2S) - 1;
      coord_1 = 0:(1/window_length_1):1;
       coord_2 = 0:(1/window_length_2):1;
       
       [wl_1,wr_1] = fdct_wrapping_window(coord_1);  
       [wl_2,wr_2] = fdct_wrapping_window(coord_2);     
       

       lowpass_1 = [wl_1, ones(1,newN1-length(wl_1)-length(wr_1)), wr_1];
       lowpass_2 = [wl_2, ones(1,newN2-length(wl_2)-length(wr_2)), wr_2];         
       lowpass = lowpass_1'*lowpass_2;
       hipass = sqrt(1 - lowpass.^2);
       Xhi(Xlow_index_1, Xlow_index_2) = Xlow .* hipass;                            
       Xlow = Xlow .* lowpass;
       
       

    % Loop: angular decomposition
    l = 0;
    nbquadrants = 2 + 2*(~is_real);              
    
    if nbangles(1,j)~=nbangles(2,j)
        notequalangular=1;
    else     notequalangular=0;
    end
    
    
    [a,b]=size(Xhi);
    vert=floor((a)/2);
    horiz=floor((b)/2);

    for quadrant = 1:nbquadrants   
            %First setting up C1 and C2 for the case of uneven angular
            %divisions between quadrants:
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
                  
            
 
            myfourMhoriz= horiz * (mod(quadrant,2)==1) + vert * (mod(quadrant,2)==0);        
            myfourMvert= vert * (mod(quadrant,2)==1) + horiz * (mod(quadrant,2)==0);        
            M_horiz=(myfourMhoriz/4);
            M_vert=(myfourMvert/4);

            nbangles_perquad = nbangles(quadrant,j);

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
            
            wedge_midpoints = (wedge_endpoints(1:(end-1)) + wedge_endpoints(2:end))/2;                   
            
            
                
        % Left corner wedge:
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
        [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
               
        for row = Y_corner
            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+1+abs(cols-1)));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols > 0);
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
            

        end;

           %wrapping_YY: Y coodinates of every single point in the
           %parallelogram representing the wedge (wrapped)
        slope_wedge_right = (floor(4*M_horiz)+1 - wedge_midpoints(1))/floor(4*M_vert);
        mid_line_right = wedge_midpoints(1) + slope_wedge_right*(wrapped_YY - 1);


        % not integers in general
        coord_right = 1/2 + floor(4*M_vert)/(wedge_endpoints(2) - wedge_endpoints(1)) * ...
            (wrapped_XX - mid_line_right)./(floor(4*M_vert)+1 - wrapped_YY);            
                        
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) = ...
            wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) + (wrapped_YY-1)/floor(4*M_vert) == 2) + 1;
  
        if notequalangular==1
           coord_corner = C1left + C2left * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        else     
            
            C2 = 1/(1/(2*(floor(4*M_horiz))/(wedge_endpoints(1) - 1) - 1)+ 1/(2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1));
    % 
            C1 = C2 / (2*(floor(4*M_vert))/(first_wedge_endpoint_vert - 1) - 1);

             coord_corner = C1left + C2left * ((wrapped_XX - 1)/(floor(4*M_horiz)) - (wrapped_YY - 1)/(floor(4*M_vert))) ./ ...
            (2-((wrapped_XX - 1)/(floor(4*M_horiz)) + (wrapped_YY - 1)/(floor(4*M_vert))));
        end
                                       
        wl_left = fdct_wrapping_window(coord_corner);
        
        [~,wr_right] = fdct_wrapping_window(coord_right);
    
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
                
                 
        
        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
                C{j}{l} = sqrt(2)*real(x);


                    C{j}{l+sum(nbangles(1,j)+nbangles(2,j))} = sqrt(2)*imag(x);

        end;
                                                                                
        
        
        
        
        %Regular wedges:    
        %This is 3*M2 (M2 is the size of the next smaller dimension: 3*M1=
        %M1+M2 (where M1 is the size of the larger dimension).

        length_wedge = floor(4*M_vert)-floor(0.5*smoothwinsize(quadrant,j-1));




                
        Y = 1:length_wedge;
        
   
        first_row = floor(4*M_vert)+2-ceil((length_wedge+1)/2)+...
            mod(length_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
            
        

        
        for subl = 2:(nbangles_perquad-1);          
            l = l+1;                       
            width_wedge = wedge_endpoints(subl+1) - wedge_endpoints(subl-1) + 1;
            slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(subl))/floor(4*M_vert) ;          
            left_line = round(wedge_endpoints(subl-1) + slope_wedge*(Y - 1));                        
            [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_wedge,width_wedge));            
            first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
                mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2)); 
            
            for row = Y
                cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);         
                new_row = 1 + mod(row - first_row, length_wedge);
                wrapped_data(new_row,:) = Xhi(row,cols);                
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
            
            wrapped_data = wrapped_data .* (wl_left .* wr_right);
            
            


            switch is_real
                case 0
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
                case 1
                    wrapped_data = rot90(wrapped_data,-(quadrant-1));
                    x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
                    C{j}{l} = sqrt(2)*real(x);

                        C{j}{l+sum(nbangles(1,j)+nbangles(2,j))} = sqrt(2)*imag(x);


                    
            end;                        
        end;
       

        

        % Right corner wedge:

        l=l+1;
        width_wedge = 4*floor(4*M_horiz) + 3 - wedge_endpoints(end) - wedge_endpoints(end-1);
        slope_wedge = ((floor(4*M_horiz)+1) - wedge_endpoints(end))/floor(4*M_vert);
        left_line = round(wedge_endpoints(end-1) + slope_wedge*(Y_corner - 1));
        first_row = floor(4*M_vert)+2-ceil((length_corner_wedge+1)/2)+...
            mod(length_corner_wedge+1,2)*(quadrant-2 == mod(quadrant-2,2));
        first_col = floor(4*M_horiz)+2-ceil((width_wedge+1)/2)+...
            mod(width_wedge+1,2)*(quadrant-3 == mod(quadrant-3,2));
         [wrapped_data, wrapped_XX, wrapped_YY] = deal(zeros(length_corner_wedge,width_wedge));

        for row = Y_corner

            cols = left_line(row) + mod((0:(width_wedge-1))-(left_line(row)-first_col),width_wedge);
            admissible_cols = round(1/2*(cols+2*floor(4*M_horiz)+1-abs(cols-(2*floor(4*M_horiz)+1))));
            new_row = 1 + mod(row - first_row, length_corner_wedge);
            wrapped_data(new_row,:) = Xhi(row,admissible_cols) .* (cols <= (2*floor(4*M_horiz)+1));
            wrapped_XX(new_row,:) = XX(row,admissible_cols);
            wrapped_YY(new_row,:) = YY(row,admissible_cols);
        end;
        
        slope_wedge_left = ((floor(4*M_horiz)+1) - wedge_midpoints(end))/floor(4*M_vert);
        mid_line_left = wedge_midpoints(end) + slope_wedge_left*(wrapped_YY - 1);
        coord_left = 1/2 + floor(4*M_vert)/(wedge_endpoints(end) - wedge_endpoints(end-1)) * ...
            (wrapped_XX - mid_line_left)./(floor(4*M_vert) + 1 - wrapped_YY);
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) = ...
        wrapped_XX((wrapped_XX - 1)/floor(4*M_horiz) == (wrapped_YY - 1)/floor(4*M_vert)) - 1;

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
        wrapped_data = wrapped_data .* (wl_left .* wr_right);
        switch is_real
            case 0
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                C{j}{l} = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
            case 1
                wrapped_data = rot90(wrapped_data,-(quadrant-1));
                x = fftshift(ifft2(ifftshift(wrapped_data)))*sqrt(numel(wrapped_data));
                C{j}{l} = sqrt(2)*real(x);
                    C{j}{l+sum(nbangles(1,j)+nbangles(2,j))} = sqrt(2)*imag(x);
        end;

        
        if quadrant < nbquadrants, Xhi = rot90(Xhi); end;
        
    end;
end;


% Coarsest wavelet level
C{1}{1} = fftshift(ifft2(ifftshift(Xlow)))*sqrt(numel(Xlow));
if is_real == 1,
    C{1}{1} = real(C{1}{1});
end;

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