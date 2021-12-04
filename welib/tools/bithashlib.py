# bits, bytes, encoding, hash
import hashlib
import struct 
import numpy as np 

# --------------------------------------------------------------------------------}
# --- Bytes / Bits / hex conversions
# --------------------------------------------------------------------------------{
# Types: bitstr, bitlist, int, int256, in256list, double, doublelist, asciistr, hex, hexstr, bytes

# ---- Bits
def int2bitstr(i,lead=False):
    # returns a string, e.g. int2bitstring(11) -> '1011'
    s="{0:b}".format(i)
    if lead:
        #return "0b"+s
        return bin(i)
    else:
        return s
def int2bitlist(i,lead=False):
    # returns a string, e.g. int2bitstring(11) -> '1011'
    s=int2bitstr(i,lead)
    return [int(x) for x in s]
def bitstr2int(s):
    return int(s,base=2)
def bitlist2int(l):
    return bitstr2int(''.join([str(x) for x in l]))

def bytes2bitstr(data):
    return ''.join( [int2bitstr(x).zfill(8) for x in list(data)])
#def bytes2bitstr(data):
#    return ''.join([str(x) for x in bytes2bitlist(data)])
def bytes2bitlist(data):
    return [int(x) for x in bytes2bitstr(data)]
#def bytes2bitlist(data):
#    """Turn the data (bytes), into a list of bits (1, 0)'s"""
#    l = len(data) * 8
#    result = [0] * l
#    pos = 0
#    for ch in data:
#        i = 7
#        while i >= 0:
#            if ch & (1 << i) != 0:
#                result[pos] = 1
#            else:
#                result[pos] = 0
#            pos += 1
#            i -= 1
#    return result

# ---- Hex
def int2hexstr(i,lead=False):
    # returns a string, e.g. int2hexstring(255) -> 'ff'
    s=hex(i)
    if not lead:
        s=s[2:]
    return s
def hex2int(sh):
    return int(sh,16)
def hex2int256list(s):
    # splitting the string
    if len(s) % 2 !=0:
        raise Exception('Please provide a hex string of length multiple of 2')
    splits=[s[i:i+2] for i in range(0, len(s), 2)]
    return [hex2int(x) for x in splits]
# --- Ascii
def asciistr2bytes(s):
    # note, the resturn will look the same since python renders the bytes the usual characters
    return s.encode('ascii') 
def int256list2asciistr(il):
    return ''.join([chr(x) for x in il])
def bytes2ascii(data):
    #return  data.decode('ascii')
    return ''.join([chr(x) for x in list(data)])
# --- ints, int256 and int256list
def bytes2int256list(B):
    return list(B)
def int256list2bytes(il):
    return bytes(il)
def int2bytes(i):
    # works for one byte only
    return struct.pack('>i',np.int32(i))
def int2562bytes(i):
    # works for one byte only
    return bytes([i])
    #return chr(i)
# --- Double
def double2bytes(d):
    # works for one byte only
    return struct.pack('d',s)
def doublelist2bytes(l):
    # works for one byte only
    if isinstance(l,np.ndarray):
        return l.tobytes()
    else:
        return b''.join([double2bytes(x) for x in l])


# --------------------------------------------------------------------------------}
# --- Base 36/ 64 
# --------------------------------------------------------------------------------{
# --- base 36
def base36encode(integer):
    chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sign = '-' if integer < 0 else ''
    integer = abs(integer)
    result=''
    while integer > 0:
        integer, remainder = divmod(integer, 36)
        result =chars[remainder]+result
    return sign+result  #appending sign to the obtained resultdef md5sum_file(fname):
def base36decode(number):
    return int(number, 36)    
# --- base 64 - Note base 64 has +/= whichi is not convenient
#import base64    
#base64.b64encode(array2md5_bytes(npv))


# --------------------------------------------------------------------------------}
# --- Hash 
# --------------------------------------------------------------------------------{
# --- md5
def array2md5_hexstr(x):
    if not isinstance(x,np.ndarray):
        x=np.ndarray(x)
    return hashlib.md5((x.tobytes())).hexdigest()
        
def array2md5_bytes(x):
    if not isinstance(x,np.ndarray):
        x=np.ndarray(x)
    return hashlib.md5((x.tobytes())).digest()
def md5sum_file(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# --- sha1
def array2sha1_hexstr(x):
    if not isinstance(x,np.ndarray):
        x=np.ndarray(x)
    return hashlib.sha1((x.tobytes())).hexdigest()
        
def array2sha1_bytes(x):
    if not isinstance(x,np.ndarray):
        x=np.ndarray(x)
    return hashlib.sha1((x.tobytes())).digest()


#+ bytes2bitstr(array2md5_bytes(npv)))
#nparray2md5_hexstring(np.array([1,2,3]))         
#nparray2md5_bytes(np.array([1,2,3]))         
##nparray2md5_md5(np.array([1,2]))[:5]
#print(hashlib.md5('hello'.encode()).hexdigest())

#print(base36encode(35))
#print(base36decode(base36encode(350)))
#import random
#f = random.random()
#
#print(''.join(map('{0:08b}'.format, map(ord, list(s)))))  # 8 bits per byte
#hash(np.array([1,2,3.555]).tostring())
#hashlib.sha1('a'.encode()).hexdigest()
#
#def h11(w):
#    return hashlib.md5(w.encode()).hexdigest()[:9]
#def h1(w):
#    return hashlib.md5(w.encode()).hexdigest()
#
#def h6(w):
#    h = hashlib.md5(w)
#    return h.digest().encode('base64')[:6]
#
#def h10(w):
#    h = hashlib.sha256(w)
#    return h.digest().encode('base64')[:6]
#
#
#def h4(w):
#    h = hashlib.md5(w)
#    return h.digest().encode('base64')[:7]    
#
#h1('aa')
#h1('aa'.encode())
#h11('aa'.encode())

#array2md5_hexstr(npv)
#len(array2md5_hexstr(npv))
#list(array2md5_bytes(npv))
#bitstr2int(bytes2bitstr(array2md5_bytes(npv)))
#base36encode(bitstr2int(bytes2bitstr(array2md5_bytes(npv))))
#len(base36encode(bitstr2int(bytes2bitstr(array2md5_bytes(npv)))))
#
#
#hashlib.sha1(npv.tobytes()).hexdigest()


##bytes2int256list(int2bytes(100))
#int2bytes(2055)
#int2bitstr(2**31-1)
#bytes2bitstr(int2bytes(2**31-1))
##bytes2bitlist(int2bytes(2**31-1))
#bytes2bitlist(int2562bytes(255))
#bitlist2int(int2bitlist(2000))
#bitstr2int(int2bitstr(2000))
#int2hex_string(97)
#double2bytes(5.12345)+double2bytes(5.123435)
#doublelist2bytes([5.12345,5.12345])
#doublelist2bytes(np.array([5.12345,5.12345]))
#doublelist2bytes(npv)
#list(int2bytes(1))
#int2bytes(10000000)
#list(int2562bytes(97))
#hex2int256list(int2hex_string(97))

#bytes2ascii(b'\xe0')
#bytes2int256list(b'\xe0')
#asciistr2bytes(int256list2asciistr([97,10,20,97]))
#int256list2bytes([97,10,20,97])
#bytes([97,10,20,97])
#chr(224)
#hex2int('e0')
#bytes(hex2int256list('e0e1e2'))
#bytes2bit_intlist(b0)
#bytes2bit_intlist(b'\x01')
#bytes2bit_intlist(b'\x10')
#bytes2bit_intlist(b'\xffff')
#bytes2bit_intlist(b'\xff\xff')
#bytes2bit_intlist(b'A')
#bytes2bit_intlist('A'.encode())
#bytes2bit_string(b'A')
#list('bca'.encode())
#int2bit_string(98)+' '+int2bit_string(99)+' '+int2bit_string(97)
#print(struct.pack('>i',2**31-1))
#print(struct.pack('>i',0))
#print(struct.pack('>i',1))
#print(struct.pack('f',5.1))
#bytearray(struct.pack('>i',2**31-1))[0]
#bytes2intlist(struct.pack('>i',2**31-1))
#list(b'@\xa333')
#list(b'\xa333')
#list(b'\xffff')
#s = 5.12345
#type(s)
#nps=np.array([s])
#npv=np.array([s,s])
#type(nps[0])
#type(npv)
#struct.pack('f',nps)
#struct.pack('f',s)
#struct.pack('f',float(s))
#list(struct.pack('f',float(s)))
#struct.pack('d',s)
#struct.pack('d',nps)
#nps.tobytes()
#list(nps.tobytes())
#len(nps.tobytes())
#struct.pack('dd',npv[0],npv[1])
#npv.tobytes()
#list(npv.tobytes())
#
#output_file = open('mybinfile', 'wb')
#npv.tofile(output_file)
#output_file.close()
#with open('mybinfile', 'rb') as f:
#    B=f.read()
#type(B)    
#B
#len(B)
#list('ab'.encode('ascii'))
#bytearray('ab'.encode('ascii'))
#print(int2bit_string(11))
#print(int2bit_string(11,lead=True))
#print(bin(11))
#print(bin(255)) # string
#print(hex(255)) # string
#print(int2hex_string(255)) # string
#print(int2hex_string(255,lead=True)) # string
#import struct
#print(struct.pack('>i',2**31-1))
#print(struct.pack('>i',0))
#print(struct.pack('>i',1))
#aa=struct.pack('>d',5.1)
#print(aa)
#type(aa)
#print(2147483647)
#print(2**31-1)
