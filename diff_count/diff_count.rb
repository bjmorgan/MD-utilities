#!/usr/bin/ruby -w

# November 23, 2011

# diff_count
# reads in tetrahedral occupations from "fort.40", and counts diffusion events
# reads in positions from "poscart.out" and cell lengths from the end of "restart.dat" so that distance travelled per event can be calculated.

class Site
  
  attr_accessor :id
  
  def self.ntet=( nlattice_sites )
    @@nlattice_sites = nlattice_sites
  end
  
  def initialize( site_number )
    @id = site_number
  end
  
  def is_lattice_site?
    if ( ( @id > 0 ) && ( @id <= @@nlattice_sites ) )
      return true
    else
      return false
    end
  end
  
  def is_interstitial?
    !is_lattice_site?
  end
  
  def copy
    return Site.new( id )
  end
  
  def is_oct?
    if ( @id == 0 )
      return true
    else
      return false
    end
  end
  
  def is_tet?
    if ( @id > @@nlatticesites )
      return true
    else
      return false
    end
  end
  
end
    
class Ion
  
  attr_accessor :new_site, :old_site, :previous_site, :step_departed, :first_oct, :first_tet
  attr_reader :id
  
  @@nions = 0
  
  def initialize( init_site )
    @@nions += 1
    @id = @@nions
    @new_site = Site.new( init_site )
    @previous_site = @new_site.copy
    @step_departed = nil
    clear_history
    update
  end
  
  def to_s
    puts "#{id}:{#{old_site.id},#{new_site.id}}"
  end
    
  def update
    @old_site = @new_site.copy
  end
    
  def clear_history
    @first_oct = false
    @first_tet = false
  end
  
end

nions = ARGV[0].to_i
ntet = ARGV[1].to_i
unless ARGV[2].nil?
  nskip = ARGV[2].to_i 
else
  nskip = 1
end

ions = Array.new
Site.ntet = ntet
ndiff_oct = 0
ndiff_tet = 0
noct = 0
noct_longt = 0
ntet = 0
ntet_longt = 0
nFrenkel_oct = 0
nFrenkel_tet = 0
nanions = ntet / 2

begin
  file = File.new("fort.40", "r")
rescue
  abort "fort.40 not found"
end

fileout = File.new("diff_count.out", "w")

step = 1
line = file.gets.split
line.shift # remove step number
line.each { |value| ions << Ion.new( value.to_i ) }
while ( tline = file.gets )
  line = tline.split
  step += 1
  fileout.puts "New Step: #{step}"
  line.shift
  line.each_with_index do |value, index| 
    ions[index].new_site.id = value.to_i
  end
  if step % nskip == 0
    ions.each do |ion|
      if ion.new_site.id != ion.previous_site.id # ion has moved to another site
        if ion.previous_site.is_lattice_site? # Frenkel pair created
          fileout.puts "Step #{step}: Ion #{ion.id} moved from lattice site #{ion.old_site.id}" 
          ion.step_departed = step 
          if ion.new_site.is_oct?
            nFrenkel_oct += 1
            ion.first_oct = true
          else
            nFrenkel_tet += 1
            ion.first_tet = true
          end
        elsif ion.new_site.is_lattice_site? # Frenkel pair annihilated
          if ion.new_site.id == ion.old_site.id # Ion has returned to original site. No net diffusion.
            fileout.puts "Step #{ion.step_departed}=>#{step}: Ion #{ion.id} returned to lattice site #{ion.old_site.id}"
          elsif ion.old_site.is_lattice_site? # Ion has moved to a new lattice site. Complete diffusion event.
            fileout.puts "Step #{ion.step_departed}=>#{step}: Ion #{ion.id} moved from lattice site #{ion.old_site.id} to #{ion.new_site.id}"
            ion.update
            ndiff_oct += 1 if ion.first_oct == true
            ndiff_tet += 1 if ion.first_tet == true
          end
          ion.clear_history
        end
      else # ion has stayed at the same site
        if ion.new_site.is_interstitial? 
          if ion.new_site.is_oct?
            noct_longt += 1
          else
            ntet_longt += 1
          end
        end
      end
      if ion.new_site.is_interstitial? # interstitial exists
        if ion.new_site.is_oct? 
          noct += 1
        else
          ntet += 1
        end
      end
      ion.previous_site = ion.new_site.copy
    end
  end
end
printf "%2.4f\t\tFrenkel pairs formed with Ag_i(Oh)\n" % [nFrenkel_oct.to_f / (step/nskip).to_f]
puts "%2.4f\t\tFrenkel pairs formed with Ag_i(Td)\n" % [nFrenkel_tet.to_f / (step/nskip).to_f]
puts "#{noct.to_f / step.to_f}\t\tn(Ag_i(Oh))"
puts "#{noct_longt.to_f / step.to_f}\t\tn(Ag_i(Oh)) sequential frames"
puts "#{ntet.to_f / step.to_f}\t\tn(Ag_i(Td))"
puts "#{ntet_longt.to_f / step.to_f}\t\tn(Ag_i(Td)) sequential frames"
printf "%.5f => %.5f pct\tOh diffusion events\n" % [ndiff_oct.to_f / (step/nskip).to_f, 100*ndiff_oct.to_f/nFrenkel_oct.to_f ]
printf "%.5f => %.5f pct\tTd diffusion events\n" % [ndiff_tet.to_f / (step/nskip).to_f , (100*ndiff_tet.to_f/nFrenkel_tet.to_f) ]
printf "%.5f => %.5f pct\tall diffusion events\n" % [( ndiff_oct + ndiff_tet ).to_f / (step/nskip).to_f , (100*( ndiff_oct + ndiff_tet).to_f / (nFrenkel_tet + nFrenkel_oct).to_f) ]
