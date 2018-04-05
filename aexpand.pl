#
#   Expand axial-turn metric solutions
#
#   U F R D B L
my @moves = ('U','F','R','D','B','L', 'A','C','E','G','H','N','O','P','M') ;
for ($i=0; $i<@moves; $i++) {
   $mind{$moves[$i]} = $i ;
}
$mind{'I'} = $mind{'N'} ;
$mind{'J'} = $mind{'O'} ;
$mind{'K'} = $mind{'P'} ;
sub expand {
   my $mov = shift ;
   my $fac = substr($mov, 0, 1) ;
   my $cnt = substr($mov, 1) ;
   return $mov if $mind{$fac} < 6 ;
   my $m1 = $moves[$mind{$fac} % 3] . $cnt ;
   my $m2 = $moves[$mind{$fac} % 3 + 3] . (1+int(($mind{$fac}-6)/3)) ;
   return $m1 . $m2 ;
}
while (<>) {
   s/([A-Z][123])/expand($1)/ge ;
   print ;
}
