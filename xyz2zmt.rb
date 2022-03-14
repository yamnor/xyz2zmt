#!/usr/bin/ruby

require 'optparse'

gazfile = nil
fhzfile = nil
fchfile = nil
zmtfile = nil
xyzfile = nil
gazmake = nil

OptionParser.new do |opt|
  opt.on('--igaz VALUE') do |v|
    gazfile = v
  end
  opt.on('--ifhz VALUE') do |v|
    fhzfile = v
  end
  opt.on('--ifch VALUE') do |v|
    fchfile = v
  end
  opt.on('--ixyz VALUE') do |v|
    xyzfile = v
  end
  opt.on('--ozmt VALUE') do |v|
    zmtfile = v
  end
  opt.on('--ogaz VALUE') do |v|
    gazfile = v
    gazmake = true
  end
  opt.parse!(ARGV)
end

BohrToAngstrom = 0.52917720859
RadianToDegree = 180.0 / Math::PI

def dot(a, b)
  w = 0.0
  for x in 0..2
    w += a[x] * b[x]
  end
  return w
end

def cross(a, b)
  v = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
  w = 0.0
  for x in 0..2
    w += v[x]**2
  end
  w = Math.sqrt(w)
  return v, w
end

def vector(a, b)
  v = Array.new(3)
  w = 0.0
  for x in 0..2
    v[x] = a[x] - b[x]
    w += v[x]**2
  end
  w = Math.sqrt(w)
  for x in 0..2
    v[x] /= w
  end
  return v, w
end

def read_fch(fch)
  geom = Array.new
  natoms = 0
  File.open(fch, "r") do |f|
    while tmp = f.gets
      if tmp =~ /^Number of atoms/
        natoms = (tmp.split)[4].to_i
      end
      if tmp =~ /^Current cartesian coordinates/
        while geom.size < natoms*3
          dat = f.gets.split
          dat.each do |x|
            geom << x.to_f
          end
        end
      end
    end
  end
  mol = Array.new(natoms)
  for n in 0..natoms-1
    mol[n] = Array.new(3)
    for i in 0..2
      mol[n][i] = geom[n*3+i]
    end
  end
  return mol
end

def read_xyz(xyz)
  geom = Array.new
  natoms = 0
  File.open(xyz, "r") do |f|
    natoms = f.gets.to_i
    f.gets
    natoms.times do
      dat = f.gets.split
      for i in 1..3
        geom << dat[i].to_f
      end
    end
  end
  mol = Array.new(natoms)
  for n in 0..natoms-1
    mol[n] = Array.new(3)
    for i in 0..2
      mol[n][i] = geom[n*3+i]
    end
  end
  return mol
end

def dist(geom, i, j)
  eij, rij = vector(geom[i], geom[j])
  return rij
end

def angl(geom, i, j, k)
  eij, rij = vector(geom[i], geom[j])
  ekj, rkj = vector(geom[k], geom[j])
  costheta = dot(eij, ekj)
  if costheta > 1.0
    costheta = 1.0
  end
  if costheta < -1.0
    costheta = -1.0
  end
  return theta = Math.acos(costheta)
end

def tors(geom, i, j, k, l)
  eji, rji = vector(geom[j], geom[i])
  ekj, rkj = vector(geom[k], geom[j])
  elk, rlk = vector(geom[l], geom[k])
  eji_ekj, eji_ekj_norm = cross(eji, ekj)
  ekj_elk, ekj_elk_norm = cross(ekj, elk)
  cosphi = dot(eji_ekj, ekj_elk) / (eji_ekj_norm * ekj_elk_norm)
  if cosphi > 1.0
    cosphi =  1.0
  end
  if cosphi < -1.0
    cosphi = -1.0
  end
  phi = Math.acos(cosphi)
  v_s, v_s_norm = cross(eji_ekj, ekj_elk)
  if v_s_norm != 0.0
    sgn = dot(v_s, ekj) / (v_s_norm * rkj)
  else
    sgn = 1.0
  end
  if (sgn > 0.0)
    sgn = +1.0
  end
  if (sgn < 0.0)
    sgn = -1.0
  end
  return phi *= sgn
end

def make_zmat(file, zmat)
  File.open(file, "w") do |f|
    for n in 0..zmat.size-1
      case zmat[n].size
      when 1
        f.puts "%-2s"%zmat[n]
      when 3
        f.puts "%-2s%3d%4s"%zmat[n]
      when 5
        f.puts "%-2s%3d%4s%3d%4s"%zmat[n]
      when 7
        f.puts "%-2s%3d%4s%3d%4s%3d%4s"%zmat[n]
      end
    end
  end
end

def make_zval(file, geom, zmat)
  File.open(file, "w") do |f|
    for n in 0..zmat.size-1
      case zmat[n].size
      when 3
        ni = n
        nj = zmat[n][1] - 1
        line = ""
        line += "%-5s %15.8f\n"%[zmat[n][2], dist(geom, ni, nj)]
        f.puts line
      when 5
        ni = n
        nj = zmat[n][1] - 1
        nk = zmat[n][3] - 1
        line = ""
        line += "%-5s %15.8f\n"%[zmat[n][2], dist(geom, ni, nj)]
        line += "%-5s %15.8f\n"%[zmat[n][4], angl(geom, ni, nj, nk) * RadianToDegree]
        f.puts line
      when 7
        ni = n
        nj = zmat[n][1] - 1
        nk = zmat[n][3] - 1
        nl = zmat[n][5] - 1
        line = ""
        line += "%-5s %15.8f\n"%[zmat[n][2], dist(geom, ni, nj)]
        line += "%-5s %15.8f\n"%[zmat[n][4], angl(geom, ni, nj, nk) * RadianToDegree]
        line += "%-5s %15.8f\n"%[zmat[n][6], tors(geom, ni, nj, nk, nl) * RadianToDegree]
        f.puts line
      end
    end
  end
end

def read_gaz(file)
  zmat = Array.new
  File.open(file, "r") do |f|
    temp = ""
    while line = f.gets
      temp += line
    end
    zmat = temp.split("\n")
    for n in 0..zmat.size-1
      zmat[n] = zmat[n].split
      case zmat[n].size
      when 1
        zmat[n] = [zmat[n][0]]
      when 3
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
      when 5
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
        zmat[n][3] = zmat[n][3].to_i
        zmat[n][4] = "A%02d"%[n+1]
      when 7
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
        zmat[n][3] = zmat[n][3].to_i
        zmat[n][4] = "A%02d"%[n+1]
        zmat[n][5] = zmat[n][5].to_i
        zmat[n][6] = "D%02d"%[n+1]
      end
    end
  end
  return zmat
end

def read_fhz(file)
  zmat = Array.new
  File.open(file, "r") do |f|
    f.gets
    f.gets
    temp = ""
    while line = f.gets
      temp += line
    end
    zmat = temp.split("\n")
    for n in 0..zmat.size-1
      zmat[n] = zmat[n].split
      case zmat[n].size
      when 2
        zmat[n] = [zmat[n][0]]
      when 3
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
      when 5
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
        zmat[n][3] = zmat[n][3].to_i
        zmat[n][4] = "A%02d"%[n+1]
      when 7
        zmat[n][1] = zmat[n][1].to_i
        zmat[n][2] = "B%02d"%[n+1]
        zmat[n][3] = zmat[n][3].to_i
        zmat[n][4] = "A%02d"%[n+1]
        zmat[n][5] = zmat[n][5].to_i
        zmat[n][6] = "D%02d"%[n+1]
        zmat[n][7] = 0
      end
    end
  end
  return zmat
end

#if fhzfile == nil && gazfile == nil || fhzfile != nil && gazfile != nil
#  raise ArgumentError
#end

if xyzfile != nil
  geom = read_xyz(xyzfile)
end

if fhzfile != nil
  zmat = read_fhz(fhzfile)
end

if gazfile != nil && gazmake != true
  zmat = read_gaz(gazfile)
end

if fchfile != nil
  geom = read_fch(fchfile)
end

if zmtfile != nil && geom != nil
  make_zval(zmtfile, geom, zmat)
end

if gazfile != nil && gazmake == true
  make_zmat(gazfile, zmat)
end
