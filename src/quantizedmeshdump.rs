use binrw::io::{Cursor,Result,Error, ErrorKind}; // A no_std reimplementation of std::io
use binrw::{
    BinRead // trait for reading
};
use geo::Vertex3;

mod geo;


// #[wasm_bindgen]
// pub fn main() -> Result<(), JsValue> {
//         //let filepath = std::env::args().nth(1).expect("No file path provided");

//     //let path = Path::new(&filepath);
//     //let mut file = File::open(&path)?;
    
//     let data = include_bytes!("../75953.terrain");
//     let mut cursor = Cursor::new(data);
//     let qm: geo::QuantizedMesh = geo::QuantizedMesh::read_le(&mut cursor).map_err(|e| JsValue::from_str(&format!("Failed to read quantized mesh: {}", e)))?;
//     print_header(&qm);


//     Ok(())
// }

// #[wasm_bindgen]
// pub fn main() -> Result<(), JsValue> {
//     // WASM-specific code
//     let data = include_bytes!("../75953.terrain");
//     let mut cursor = Cursor::new(data);
//     let qm: QuantizedMesh = QuantizedMesh::read_le(&mut cursor).map_err(|e| JsValue::from_str(&format!("Failed to read quantized mesh: {}", e)))?;
//     print_header(&qm);

//     Ok(())
// }

fn main() -> Result<()> {
    // Native code

    use std::{fs::File, path::Path};
    let filepath = std::env::args().nth(1).expect("No file path provided");
    let path = Path::new(&filepath);
    let mut file = File::open(&path)?;
    //let mut cursor = Cursor::new(file);
    //let qm: geo::QuantizedMesh = geo::QuantizedMesh::read_le(&mut file).unwrap();

    // let qm: geo::QuantizedMesh = match geo::QuantizedMesh::read_le(&mut cursor) {
    //     Ok(mesh) => mesh,
    //     Err(e) => return Err(" {:?}", e))),
    // };

    
    
    
    let qmheader: geo::QuantizedMeshHeader = match geo::QuantizedMeshHeader::read(&mut file) {
        Ok(mesh) => mesh,
        Err(e) => return Err(Error::new(ErrorKind::Other, format!("Failed to read quantized mesh: {:?}", e)))
    };
    println!("{qmheader:#?}");


     let vertex_data: geo::VertexData = match geo::VertexData::read(&mut file) {
        Ok(mesh) => mesh,
        Err(e) => return Err(Error::new(ErrorKind::Other, format!("Failed to read quantized mesh: {:?}", e)))
     };
     let uvh: Vec<Vertex3<&i32>> = vertex_data.as_uvh();

    println!("vertex count: {:#?}", uvh.len());

    println!("index_data_16 len {:#?}", vertex_data.index_data_short.len());
    println!("index_data_32 len {:#?}", vertex_data.index_data_long.len());

     



    //map_err(|e| binrw::io::Error::new(binrw::io::ErrorKind::Other, e))?;

    Ok(())
}



