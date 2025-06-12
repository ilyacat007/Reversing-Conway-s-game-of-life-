use std::collections::HashMap;


fn get_up(x:&[bool], n:usize) -> Vec<bool>{
    return x[..n*2].to_vec();
}

fn get_down(x:&[bool], n:usize) -> Vec<bool>{
    return x[x.len()-n*2..].to_vec();
}

fn get_left(x: &[bool], n: usize) -> Vec<bool> {
    let mut result = Vec::with_capacity(2 * n);
    for i in 0..n {
        let start = i * n;
        result.push(x[start]);
        result.push(x[start + 1]);
    }
    result
}

fn get_right(x: &[bool], n: usize) -> Vec<bool> {
    let mut result = Vec::with_capacity(2 * n);
    for i in 0..n {
        let row_start = i * n;
        let last_col = row_start + n - 1;
        let penultimate_col = last_col - 1;

        result.push(x[penultimate_col]);
        result.push(x[last_col]);
    }
    result
}

fn generate_cells(center: bool, n: u32) -> Vec<[bool; 9]> {
    let min_n = if center { 1 } else { 0 };
    let max_n = if center { 9 } else { 8 };

    if n < min_n || n > max_n {
        return Vec::new();
    }


    let non_center_indices = [0, 1, 2, 3, 5, 6, 7, 8];
    let mut result = Vec::new();

    for bits in 0u8..=255 {
        if bits.count_ones() == n {
            let mut arr = [false; 9];
            arr[4] = center;
            for bit_index in 0..8 {
                if (bits >> bit_index) & 1 == 1 {
                    let pos = non_center_indices[bit_index];
                    arr[pos] = true;
                }
            }
            result.push(arr);
        }
    }

    result
}

fn alive_variants() -> Vec<[bool; 9]> {
    let mut alive = Vec::<[bool; 9]>::new();

    alive.append(&mut generate_cells(false, 3));
    alive.append(&mut generate_cells(true, 3));
    alive.append(&mut generate_cells(true, 2));

    return alive;
}



fn dead_variants() -> Vec<[bool; 9]> {
    let mut dead = Vec::<[bool; 9]>::new();

    dead.append(&mut generate_cells(true, 0));

    dead.append(&mut generate_cells(true, 1));
    dead.append(&mut generate_cells(true, 4));
    dead.append(&mut generate_cells(true, 5));
    dead.append(&mut generate_cells(true, 6));
    dead.append(&mut generate_cells(true, 7));
    dead.append(&mut generate_cells(true, 8));


    dead.append(&mut generate_cells(false, 1));
    dead.append(&mut generate_cells(false, 2));

    return dead;
}



fn merge(first:&Vec<[i32; 4]>, second:&Vec<[i32; 4]>, third:&Vec<[i32; 4]>, fourth:&Vec<[i32; 4]>, any:i32) -> Vec<[[i32; 4]; 4]>{
    let mut answer:Vec<[[i32; 4]; 4]> = Vec::<[[i32; 4]; 4]>::new();
    let mut i = 0;

    for fi in first.iter(){
        for se in second.iter() {
            if fi[1] != se[0]{
                continue
            }

            for th in third.iter(){
                if fi[3] != th[2]{
                    continue
                }

                for fo in fourth.iter(){
                    if fo[0] != th[1]{
                        continue
                    }

                    if fo[2] != se[3] {
                        continue
                    }


                    answer.push([fi.clone(), se.clone(), th.clone(), fo.clone()]);

                    if any != 0 {
                        i += 1;
                        if i>=any {
                            return answer;
                        }
                    }

                }
            }
        }
    }

    return answer

}
fn hash_variants(types:Vec<Vec<[bool;9]>>) -> (Vec<Vec<[i32; 4]>>, HashMap<[i32; 4], Vec<bool>>) {
    let mut map:HashMap<Vec<bool>, i32> = HashMap::new();
    let mut reverse_map:HashMap<i32, Vec<bool>> = HashMap::new();
    let mut decode_map:HashMap<[i32; 4], Vec<bool>> = HashMap::new();
    
    let mut answer = vec![];

    let mut i = 0;

    for t in types{
        let mut thashes = Vec::<[i32; 4]>::new();
        
        for var in t.iter(){
            let left = get_left(var, 3);
            let right = get_right(var, 3);
            let up = get_up(var, 3);
            let down = get_down(var, 3);
            
            let l = map.entry(left.clone()).or_insert(i).clone();
            reverse_map.entry(l).or_insert(left);
            i+=1;

            let r = map.entry(right.clone()).or_insert(i).clone();
            reverse_map.entry(r).or_insert(right);
            i+=1;

            let u = map.entry(up.clone()).or_insert(i).clone();
            reverse_map.entry(u).or_insert(up);
            i+=1;

            let d = map.entry(down.clone()).or_insert(i).clone();
            reverse_map.entry(d).or_insert(down);
            i+=1;
            
            decode_map.entry([l, r, u, d]).or_insert(var.to_vec());
            thashes.push([l, r, u, d]);
        }
        answer.push(thashes)
    }
    
    return (answer, decode_map)
}

fn hash_variants2(types:Vec<Vec<[[i32;4]; 4]>>, map: &mut HashMap<[i32;2], i32>,  decode_map: &mut HashMap<[i32; 4], [[i32; 4];4]>) -> Vec<Vec<[i32; 4]>> {

    let mut answer = vec![];
    let mut i = 0;

    for t in types{
        let mut thashes = Vec::<[i32; 4]>::new();

        for var in t.iter(){
            let left = [var[0][0], var[2][0]];
            let right = [var[1][1], var[3][1]];
            let up = [var[0][2], var[1][2]];
            let down = [var[2][3], var[3][3]];

            let l = map.entry(left).or_insert(i).clone();
            i+=1;

            let r = map.entry(right).or_insert(i).clone();
            i+=1;

            let u = map.entry(up).or_insert(i).clone();
            i+=1;

            let d = map.entry(down).or_insert(i).clone();
            i+=1;

            decode_map.entry([l, r, u, d]).or_insert(var.clone());
            thashes.push([l, r, u, d]);
        }
        answer.push(thashes)
    }

    return answer
}




fn illustrate(x:&Vec<bool>){
    let n = x.len().isqrt();
    for i in x.chunks(n) {
        for j in i{
            if *j{
                print!("1 ")
            }else {
                print!("0 ")
            }
        }
        println!();
    }
    println!();
}





fn reconstruct(data: [[i32;4];4], l:usize, init_map:&HashMap<[i32; 4], Vec<bool>>, leveled_maps: &mut HashMap<usize, HashMap<[i32; 4], [[i32; 4];4]>>) -> Vec<bool> {
    let mut answ = vec![false; (l+2).pow(2)];

    for n in 0..data.len(){
        let tile = match l {
            2 => {init_map.get(&data[n]).unwrap().clone()},
            _ => {
                let level = leveled_maps.get(&(l/2)).unwrap();
                let d= level.get(&data[n]).unwrap().clone();
                reconstruct(d, l/2, init_map, leveled_maps)
            }
        };


        let i0 = (n/2) * ((l/2) * (l + 2)) + (n%2)*(l/2);

        for i in 0..tile.len(){
            let ind = (l/2)*(i/(l/2 + 2)) + i0 + i;
            answ[ind] = tile[i];
        }
    }

    return answ;


}



fn step(current: &Vec<bool>) -> Vec<bool> {
    let n = current.len().isqrt();
    let mut next = vec![false; n * n];

    for i in 0..n {
        for j in 0..n {
            let idx = i * n + j;
            let mut live_neighbors = 0;

            for di in [-1, 0, 1] {
                for dj in [-1, 0, 1] {
                    if di == 0 && dj == 0 { continue; }

                    let ni = i as isize + di;
                    let nj = j as isize + dj;

                    if ni >= 0 && ni < n as isize && nj >= 0 && nj < n as isize {
                        let nidx = ni as usize * n + nj as usize;
                        if current[nidx] {
                            live_neighbors += 1;
                        }
                    }
                }
            }

            next[idx] = match (current[idx], live_neighbors) {
                (true, 2) | (true, 3) => true,
                (false, 3) => true,
                _ => false,
            };
        }
    }

    next
}

fn split_left_right(input:Vec<bool>, l:usize) -> (Vec<bool>, Vec<bool>) {
    let mut left = vec![];
    let mut right = vec![];

    let mut i = 0;

    for ch in input.chunks(l)  {
        if i%2==0 {
            left.append(&mut ch.to_vec());
        }else {
            right.append(&mut ch.to_vec())
        }
        i+=1
    }

    return (left, right)
}

fn split_up_down(input:Vec<bool>, l:usize) -> (usize, Vec<bool>, Vec<bool>){
    let l2 = l/2;

    let up = input[..l*l2].to_vec();
    let down = input[l*l2..].to_vec();

    return (l2, up, down)
}

fn go(input:Vec<bool>, dead: &Vec<[i32; 4]>, alive:&Vec<[i32; 4]>, forward_maps: HashMap<usize, HashMap<[i32; 2], i32>>, backward_maps: HashMap<usize, HashMap<[i32; 4], [[i32; 4]; 4]>>) ->( Vec<[[i32; 4]; 4]>, HashMap<usize, HashMap<[i32; 2], i32>>, HashMap<usize, HashMap<[i32; 4], [[i32; 4]; 4]>>){
    let l = input.len().isqrt();

    let (l2, up, down) = split_up_down(input, l);

    let (left_up, right_up) = split_left_right(up, l2);
    let (left_down, right_down) = split_left_right(down, l2);


    if l2 == 1{
        let mut ch = vec![];
        for i in [left_up[0], right_up[0], left_down[0], right_down[0]] {
            if i {
                ch.push(alive)
            }
            else {
                ch.push(dead)
            }
        }

        return (merge(&ch[0], &ch[1], &ch[2], &ch[3], 0), forward_maps, backward_maps)

    }else {
        let (lu, forward_maps, backward_maps)  = go(left_up, dead, alive, forward_maps, backward_maps);
        let (ru, forward_maps, backward_maps) = go(right_up, dead, alive, forward_maps, backward_maps);
        let (ld, forward_maps, backward_maps) = go(left_down, dead, alive, forward_maps, backward_maps);
        let (rd, mut forward_maps, mut backward_maps) = go(right_down, dead, alive, forward_maps, backward_maps);

        let map:HashMap<[i32;2], i32> = HashMap::new();
        let decode_map:HashMap<[i32; 4], [[i32; 4];4]> = HashMap::new();

        let forw= forward_maps.entry(l2).or_insert(map);
        let backw= backward_maps.entry(l2).or_insert(decode_map);


        println!("solving size: {l}");
        let types_variants = hash_variants2( vec![lu, ru, ld, rd], forw, backw);



        let mut any = 0;
        if l == 8 {
            any = 1;
        }

        return (merge(&types_variants[0], &types_variants[1], &types_variants[2], &types_variants[3], any), forward_maps, backward_maps)
    }

}

fn main() {
    let dead = dead_variants();
    let alive = alive_variants();

    let (types_variants, decode0_map) = hash_variants(vec![dead, alive]);


    // let input = vec![true, false, false, false,
    //                             false, true, false, false,
    //                             false, true, false, false,
    //                             false, false, false, false,];

    let input = vec![false, false, false, false, false, false, false, false,
                               false, false, false, false, false, false, false, false,
                               false, false, false, true, true, false, false, false,
                               false, false, true, true, true, false, false, false,
                               false, true, false, true, false, true, false, false,
                               false, false, true, false, true, false, false, false,
                               false, false, true, false, true, false, false, false,
                               false, false, true, false, true, false, false, false,];

    let l = input.len().isqrt();


    let forward_maps:HashMap<usize, HashMap<[i32;2], i32>> = HashMap::new();
    let backward_maps:HashMap<usize, HashMap<[i32; 4], [[i32; 4];4]>> = HashMap::new();

    let (answer, _, mut backward_maps) = go(input.clone(), &types_variants[0], &types_variants[1], forward_maps, backward_maps);//merge(&types_variants[0], &types_variants[1], &types_variants[2], &types_variants[3]);

    println!("backward maps: {:?}", backward_maps.len());
    println!("found decisions: {:?}", answer.len());

    let res = reconstruct(answer[0], l, &decode0_map, &mut backward_maps);

    illustrate(&res);
    illustrate(&step(&res));
    illustrate(&input);

}
