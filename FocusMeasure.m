function Image_FM = FocusMeasure(Image, Method)
%***** this function computes the different focus measure used to
%  detect the sharpest focus for DFF application in Light field
%  the following focus measure are currently implemented
%
%  1: SML: sum of modified laplacian
%  2: IME: image energy
%  3: RDF
% h_f= fspecial("gaussian",10,2);
switch Method
    case 'SML'
        hx = [-1,2,-1]/4;
        Image_xx = conv2(hx, 1, Image,'same');
        Image_yy = conv2(1, hx, Image,'same');
        Image_FM = abs(Image_xx) + abs(Image_yy);
%         h = ones(2,2)./4;
%         Image_IME = Image.^2;
%         Image_IME = conv2(Image_IME, h,'same');
%         Image_FM = Image_FM .* Image_IME.^2;
        
%         Image_FM = imfilter(Image_FM,h_f,'same');
%        Image_FM = imgaussfilt(Image_FM,2);
    case 'RSML'
        hx = [-1,0,2,0,-1]/4;
        Image_xx = conv2(hx, 1, Image,'same');
        Image_yy = conv2(1, hx, Image,'same');
        Image_FM = abs(Image_xx) + abs(Image_yy);
%         h = ones(3,3)./9;
%         Image_IME = Image.^2;
%         Image_IME = conv2(Image_IME, h,'same');
%         Image_FM = Image_FM .* Image_IME.^2;
%         Image_FM = imgaussfilt(Image_FM,3);

    case 'RSML2'
        hx = [-1,-1,0,4,0,-1,-1]/2;
        Image_xx = conv2(hx, 1, Image,'same');
        Image_yy = conv2(1, hx, Image,'same');
        Image_FM = abs(Image_xx) + abs(Image_yy);
%         Image_FM = imgaussfilt(Image_FM,3);

    case 'RSML3'
        hx = [-1,0, -1,0,-2,0,8,0,-2,0,-1,0,-1]/8;
        Image_xx = conv2(hx, 1, Image,'same');
        Image_yy = conv2(1, hx, Image,'same');
        Image_FM = abs(Image_xx) + abs(Image_yy);
%         Image_FM = imgaussfilt(Image_FM,3);

    case 'IME'
        h = ones(5,5)./25;
        Image_IME = Image.^2;
        Image_FM = conv2(Image_IME, h,'same');

    case 'SPARC'
        fun = @(x) nnz(x(:));
        Image_FM = -nlfilter(Image, [15,15],fun);

    case 'RDF'
        h = zeros(15,15);
        x = (1:15)-8;
        y = x;
        [x,y] = meshgrid(x,y);
        inner_ring = ((x.^2+y.^2)<=2^2);
        mask_ring =  ((x.^2+y.^2)>4^2) & (((x.^2+y.^2)<=8^2));
        h(mask_ring) = -2./nnz(mask_ring);
        h(inner_ring) = 2./nnz(inner_ring);
        Image_FM = conv2(Image,h,'same');
%         Image_FM = imgaussfilt(Image_FM,5);

    otherwise
        h = ones(5,5);
        Image_IME = Image.^2;
        Image_FM = conv2(Image_IME, h,'same');
end

end
