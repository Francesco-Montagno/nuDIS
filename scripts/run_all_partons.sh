#!/bin/bash

template_card="card/events_card.dat"

processes=("d_p" "s_p" "ubar_p" "cbar_p" )

for p in "${processes[@]}"; do
  output_card="card/run_${p}.dat"
  echo "Creating run card $output_card with process = $p"

  # Usa awk per sostituire la riga process= con il processo corrente
  awk -v proc="$p" '
    /^process[[:space:]]*=/ {
      print "process       = " proc
      next
    }
    { print }
  ' "$template_card" > "$output_card"

  echo "Running python3 ./events.py with run card $output_card"
  python3 ./events.py -c "$output_card"

  echo "Removing temporary run card $output_card"
  rm "$output_card"
done

echo "All processes completed."
