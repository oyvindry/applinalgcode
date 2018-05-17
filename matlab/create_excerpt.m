function img=create_excerpt()
    img = double(imread('images/lena.png','png'));
    img = img(129:384,129:384,:);
