BEGIN{FS=OFS="\t"} {
    if ($2==0) {
        $2=99;
        if (($4+$9-r) > $4) $8=$4+$9-r; else $8=$4;
    } else {
        $2=147;
        if (($4-$9+r) < $4) $8=$4-$9+r; else $8=$4;
        $9=$9*-1;
    }
    $7="="; $11=$10;
} {gsub(/./, "F", $11)} {print}