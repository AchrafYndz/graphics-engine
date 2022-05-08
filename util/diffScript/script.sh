for f in $1/*.ini
  do
    y=${f%.ini}
    p="$y.png"
    b="$y.bmp"
    echo "Comparing $p and $b"
    ../../cmake-build-debug/engine $f

    compare -verbose -metric AE $p $b temp.png

    read -p "Press any key to continue" mainmenuinput
done
