using Images;
h = 5;

dir = ARGS[1][length(ARGS[1])] == '/' ? ARGS[1] : ARGS[1] * "/";
files = readdir(dir);

for file in files
  if contains(file, ".jpg")
    img = Gray.(load(dir*file));
    (max_height, max_width) = size(img);

    template = ones(h, max_width);
    vs = map(i->norm(template - img[i:i+h-1, :]), 1:max_height-h);

    inds = find(x->x<0.1, vs[1:Int(floor(max_height/3))]);
    if length(inds) != 0
      ind = inds[length(inds)];
      out = img[ind:end, :];
      save(dir*file, out);
    else
      rm(dir*file)
    end
  end
end